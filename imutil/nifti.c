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

/* Helper function to read a NIFTI image (.nii, .nii.gz).
 * Prior to calling this function, use init_im(im).
 * This function allocates memory.
 */
int read_nii(const char *path, Image *const im)
{

	nifti_image *nifti;
	int x, y, z, i, dim_counter;

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

        // Check the dimensionality
	if (dim_counter > 3) {
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
	im->nc = 1;
	im_default_stride(im);
	im_resize(im);

#define IM_COPY_FROM_TYPE(type) \
    SIFT3D_IM_LOOP_START(im, x, y, z)   \
        SIFT3D_IM_GET_VOX(im, x, y, z, 0) = (float) ( \
                (double) ((type *)nifti->data)[\
                        SIFT3D_IM_GET_IDX(im, x, y, z, 0)] * \
                (double) nifti->scl_slope + \
                (double) nifti->scl_inter); \
    SIFT3D_IM_LOOP_END

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
	size_t i;

	const int dims[] = { 3, im->nx, im->ny, im->nz, 0, 0, 0, 0 };

	// Verify inputs
	if (im->nc != 1) {
		SIFT3D_ERR("write_nii: unsupported number of "
			"channels: %d. This function only supports single-"
			"channel images.", im->nc);
		return SIFT3D_FAILURE;
	}

	// Init a nifti struct and allocate memory
	if ((nifti = nifti_make_new_nim(dims, DT_FLOAT32, 1))
	    == NULL)
		goto write_nii_quit;

        // Copy the units
        nifti->dx = im->ux;
        nifti->dy = im->uy;
        nifti->dz = im->uz;

	// Copy the data
	for (i = 0; i < im->size; i++) {
		((float *)nifti->data)[i] = im->data[i];
	}

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
