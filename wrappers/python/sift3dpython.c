/*
 *  C library routines for the Python wrapper
 */

#include <stdlib.h>
#include <string.h>

#include "imutil.h"
#include "immacros.h"

/* A simpler version of the Image struct, for Python compatibility */
typedef struct _Python_Image {
	float *data;		// Raster of voxel values ~16MB
	int dims[3];		// Dimensions in x, y, and z
        double units[3];        // The units
        int nc;                 // The number of channels
        int error_code;         // Success of the last operation
} Python_Image;

/* Internal functions */
static void python_init_im(Python_Image *const py);
static int im2py(const Image *const im, Python_Image *const py);
static int verify_default_stride(const Image *const im);

/* Initialize a Python_Image before it can be used. */
static void python_init_im(Python_Image *const py) {

        int i;

        py->data = NULL; 
        py->error_code = SIFT3D_SUCCESS; 
}

/* Main function to read an image. Returns a struct which must be freed later. */
Python_Image python_im_read(const char *path) {

        Python_Image py;
        Image im;

        // Initialize the images
        init_im(&im);
        python_init_im(&py);

        // Read the image
        if (py.error_code = im_read(path, &im))
                goto im_read_quit;

        // Shallow copy the data
        if (py.error_code = im2py(&im, &py))
                goto im_read_quit;

        return py;

im_read_quit:
        im_free(&im);
        return py;
}

/* Frees the data in a Python_Image. */
void python_im_free(Python_Image *const py) {

        /* Free the memory */
        if (py != NULL && py->data != NULL) {
                free(py->data);
                py->data = NULL;
        }

        /* Mark the operation as a success */
        py->error_code = SIFT3D_SUCCESS;
}

/* Wrapper for the native function */
const char *python_im_read_get_error(const int code) {
        im_read_get_error(code);
}

/* Verifies that an image has the default strides. */
static int verify_default_stride(const Image *const im) {

        Image dummy;
        int ret;

        /* Make a dummy image to get the default strides */
        init_im(&dummy);
	memcpy(SIFT3D_IM_GET_DIMS(&dummy), SIFT3D_IM_GET_DIMS(im), 
                IM_NDIMS * sizeof(int));
        dummy.nc = im->nc;
        im_default_stride(&dummy);

        /* Test that the image uses the default strides . */
        ret = memcmp(SIFT3D_IM_GET_STRIDES(im), SIFT3D_IM_GET_STRIDES(&dummy), 
                IM_NDIMS * sizeof(size_t));

        /* Free the dummy (Does nothing) */
        im_free(&dummy);

        return ret;
}

/* Shallow copy an Image struct into a simpler, Python-friendly version. Does
 * not copy the image data. Image must be in the default stride. */
static int im2py(const Image *const im, Python_Image *const py) {

        int i;

        /* Initi

        /* Verify the strides */
        if (verify_default_stride(im))
                return SIFT3D_FAILURE;

        /* Shallow copy the fields */
        py->data = im->data;
        memcpy(py->dims, SIFT3D_IM_GET_DIMS(im), IM_NDIMS * sizeof(int));
        memcpy(py->units, SIFT3D_IM_GET_UNITS(im), IM_NDIMS * sizeof(double));
        py->nc = im->nc;

        return SIFT3D_SUCCESS;
}
