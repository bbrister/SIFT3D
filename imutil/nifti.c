/* -----------------------------------------------------------------------------
 * nifti.c 
 * -----------------------------------------------------------------------------
 * Copyright (c) 2017 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Implementation of the nifticlib wrapper for reading and writing NIFTI images.
 * -----------------------------------------------------------------------------
 */

/* SIFT3D includes */ 
#include "imutil.h"
#include "immacros.h"
#include "nifti.h"

#ifndef SIFT3D_WITH_NIFTI
/* Return error messages if this was not compiled with NIFTI support. */ 

static int nii_error_message() {
        SIFT3D_ERR("nii_error_message: SIFT3D was not compiled with NIFTI "
                "support!\n");
        return SIFT3D_WRAPPER_NOT_COMPILED;
}

int read_nii(const char *path, Image *const im) {
        return nii_error_message();
}

int write_nii(const char *path, const Image *const im) {
        return nii_error_message();
}

#else

/* Standard includes */
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/* Nifti includes */
#include <nifti1_io.h>

/* Macro returning NIFTI data index for the channel dimension */
#define SIFT3D_NIM_GET_IDX(im, x, y, z, c) ( \
        (x) + (y) * (im)->nx + (z) * (im)->nx * (im)->ny + \
        (c) * (im)->nx * (im)->ny * (im)->nz )

/* Helper function to read a NIFTI image (.nii, .nii.gz).
 * Prior to calling this function, use init_im(im).
 * This function allocates memory.
 */
int read_nii(const char *path, Image *const im)
{

	nifti_image *nifti;
        double slope;
	int x, y, z, c, i, dim_counter;

        const size_t nim_channel_stride = im->nx * im->ny * im->nz;

	// Read NIFTI file
	if ((nifti = nifti_image_read(path, 1)) == NULL) {
		SIFT3D_ERR("read_nii: failure loading file %s", path);
                return SIFT3D_FAILURE;
	}

	// Find the dimensionality of the array, given by the last dimension
	// greater than 1. Note that the dimensions begin at dim[1].
	for (dim_counter = nifti->ndim; dim_counter > 0; dim_counter--) {
		if (nifti->dim[dim_counter] > 1) {
			break;
		}
	}

        // Check the dimensionality. 4D is interpreted as a 3D array with
        // multiple channels.
	if (dim_counter > 4) {
		SIFT3D_ERR("read_nii: file %s has unsupported "
			"dimensionality %d\n", path, dim_counter);
		goto read_nii_quit;
	}

        // Fill the trailing dimensions with 1
        for (i = dim_counter; i < IM_NDIMS; i++) {
                SIFT3D_IM_GET_DIMS(im)[i] = 1;
        }

	// Store the real world coordinates
	im->ux = nifti->dx;
	im->uy = nifti->dy;
	im->uz = nifti->dz;

	// Resize im    
	im->nx = nifti->nx;
	im->ny = nifti->ny;
	im->nz = nifti->nz;
	im->nc = dim_counter == 4 ? nifti->nt : 1;
	im_default_stride(im);
	im_resize(im);

        // Ignore the slope if it's zero. This is an ill-formatted image.
        slope = nifti->scl_slope;
        if (slope == 0.0) slope = 1.0;

        // Macro to copy the data for each type
#define IM_COPY_FROM_TYPE(type) \
        SIFT3D_IM_LOOP_START_C(im, x, y, z, c)   \
                SIFT3D_IM_GET_VOX(im, x, y, z, c) = (float) ( \
                        (double) ((type *) nifti->data)[\
                                SIFT3D_NIM_GET_IDX(im, x, y, z, c)] * \
                                (double) slope + (double) nifti->scl_inter); \
        SIFT3D_IM_LOOP_END_C

	// Copy the data into im, applying the slope and intercept
	switch (nifti->datatype) {
	case NIFTI_TYPE_UINT8:
		IM_COPY_FROM_TYPE(uint8_t);
		break;
	case NIFTI_TYPE_INT8:
		IM_COPY_FROM_TYPE(int8_t);
		break;
	case NIFTI_TYPE_UINT16:
		IM_COPY_FROM_TYPE(uint16_t);
		break;
	case NIFTI_TYPE_INT16:
		IM_COPY_FROM_TYPE(int16_t);
		break;
	case NIFTI_TYPE_UINT32:
		IM_COPY_FROM_TYPE(uint32_t);
		break;
	case NIFTI_TYPE_INT32:
		IM_COPY_FROM_TYPE(int32_t);
		break;
	case NIFTI_TYPE_UINT64:
		IM_COPY_FROM_TYPE(uint64_t);
		break;
	case NIFTI_TYPE_INT64:
		IM_COPY_FROM_TYPE(int64_t);
		break;
	case NIFTI_TYPE_FLOAT32:
		IM_COPY_FROM_TYPE(float);
		break;
	case NIFTI_TYPE_FLOAT64:
		IM_COPY_FROM_TYPE(double);
		break;
	case NIFTI_TYPE_FLOAT128:
	case NIFTI_TYPE_COMPLEX128:
	case NIFTI_TYPE_COMPLEX256:
	case NIFTI_TYPE_COMPLEX64:
	default:
		SIFT3D_ERR("read_nii: unsupported datatype %s \n",
			nifti_datatype_string(nifti->datatype));
                goto read_nii_quit;
	}
#undef IM_COPY_FROM_TYPE

	// Clean up NIFTI data
	nifti_free_extensions(nifti);
	nifti_image_free(nifti);

	return SIFT3D_SUCCESS;

read_nii_quit:
        nifti_free_extensions(nifti);
        nifti_image_free(nifti);
	return SIFT3D_FAILURE;
}

/* Write a Image to the specified path, in NIFTI format.
 * The path extension must be one of (.nii, .nii.gz). */
int write_nii(const char *path, const Image *const im)
{

	nifti_image *nifti;
        int x, y, z, c;

        const size_t nim_channel_stride = im->nx * im->ny * im->nz;
        const int multi_channel = im->nc > 1;
	const int dims[] = {multi_channel ? 4 : 3, 
            im->nx, im->ny, im->nz, multi_channel ? im->nc : 0, 0, 0, 0};

	// Initialize a nifti struct and allocate memory
	if ((nifti = nifti_make_new_nim(dims, DT_FLOAT32, 1)) == NULL)
		goto write_nii_quit;

        // Set the slope and intercept to do nothing
        nifti->scl_slope = 1.0;
        nifti->scl_inter = 0.0;

        // Copy the units
        nifti->dx = im->ux;
        nifti->dy = im->uy;
        nifti->dz = im->uz;
        if (multi_channel) 
            nifti->dt = 0.f; // Channels have no size

	// Copy the data
        SIFT3D_IM_LOOP_START_C(im, x, y, z, c)
                ((float *) nifti->data)[SIFT3D_NIM_GET_IDX(im, x, y, z, c)] =
                        SIFT3D_IM_GET_VOX(im, x, y, z, c);
        SIFT3D_IM_LOOP_END_C

	if (nifti_set_filenames(nifti, path, 0, 1))
		goto write_nii_quit;

	// Sanity check
	if (!nifti_nim_is_valid(nifti, 1))
		goto write_nii_quit;

	nifti_image_write(nifti);
	nifti_free_extensions(nifti);
	nifti_image_free(nifti);

	return SIFT3D_SUCCESS;

 write_nii_quit:
	if (nifti != NULL) {
		nifti_free_extensions(nifti);
		nifti_image_free(nifti);
	}
	return SIFT3D_FAILURE;
}

#endif
