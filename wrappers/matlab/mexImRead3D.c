/* -----------------------------------------------------------------------------
 * mexImRead3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2017 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex interface to read 3D images.
 * -----------------------------------------------------------------------------
 */

#include "imutil.h"
#include "mexutil.h"
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

        Image im;
        const mxArray *mxPath, *mxIm;
        const char *path;

/* Clean up and print an error */
#define CLEAN_AND_QUIT(name, msg, expected) { \
                im_free(&im); \
                if (expected) { \
                        err_msg(name, msg); \
                } else { \
                        err_msgu(name, msg); \
                } \
        }

	// Verify the number of inputs
	if (nrhs != 1)
                err_msgu("main:numInputs", "This function takes 1 input.");

        // Verify the number of outputs
        if (nlhs > 2) 
                err_msgu("main:numOutputs", "This function takes 2 outputs.");

        // Assign the inputs
        mxPath = prhs[0];

        // Initialize intermediates
        init_im(&im);

        // Get the path string
        if ((path = mxArrayToString(mxPath)) == NULL)
                CLEAN_AND_QUIT("main:getPath", "Failed to convert the input "
                        "to a string", SIFT3D_FALSE);

        // Load the image
        switch (im_read(path, &im)) {
                case SIFT3D_SUCCESS:
                        break;
                case SIFT3D_FILE_DOES_NOT_EXIST:
                        CLEAN_AND_QUIT("main:dne", "File does not exist", 
                                SIFT3D_TRUE);
                case SIFT3D_UNSUPPORTED_FILE_TYPE:
                        CLEAN_AND_QUIT("main:unsupportedType", "Unsupported file "
                                "type", SIFT3D_TRUE);
                case SIFT3D_WRAPPER_NOT_COMPILED:
                        CLEAN_AND_QUIT("main:notCompiled", "Recompile SIFT3D "
                                "with support for this file type", SIFT3D_TRUE);
                case SIFT3D_UNEVEN_SPACING:
                        CLEAN_AND_QUIT("main:unevenSpacing", "The series has "
                        "uneven slice spacing", SIFT3D_TRUE);
                case SIFT3D_INCONSISTENT_AXES:
                        CLEAN_AND_QUIT("main:inconsistentAxes", "The series "
                        "has inconsistent slice axes", SIFT3D_TRUE);
                case SIFT3D_DUPLICATE_SLICES:
                        CLEAN_AND_QUIT("main:duplicateSlices", "The series "
                        "contains slices in duplicate locations", SIFT3D_TRUE);
                default:
                        CLEAN_AND_QUIT("main:unexpected", "Unexpected error "
                                "reading the image", SIFT3D_TRUE);
        }

        // Convert the output image to a MATLAB array
        if ((plhs[0] = im2mx(&im)) == NULL)
                CLEAN_AND_QUIT("main:im2mx", "Failed to convert image to an "
                        "mxArray", SIFT3D_FALSE);

        // Convert the output units to a MATLAB array
        if ((plhs[1] = units2mx(&im)) == NULL)
                CLEAN_AND_QUIT("main:im2mx", "Failed to convert units to an "
                        "mxArray", SIFT3D_FALSE);

        // Clean up
        im_free(&im);

#undef CLEAN_AND_QUIT
}

