/* -----------------------------------------------------------------------------
 * mexImWrite3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex interface to write 3D images.
 * -----------------------------------------------------------------------------
 */

#include "imutil.h"
#include "mexutil.h"
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

        Image im;
        const mxArray *mxPath, *mxIm, *mxUnits;
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
	if (nrhs != 3)
                err_msgu("main:numInputs", "This function takes 3 inputs.");

        // Verify the number of outputs
        if (nlhs != 0) 
                err_msgu("main:numOutputs", "This function takes no outputs.");

        // Assign the inputs
        mxPath = prhs[0];
        mxIm = prhs[1];
        mxUnits = prhs[2];

        // Initialize intermediates
        init_im(&im);

        // Get the path string
        if ((path = mxArrayToString(mxPath)) == NULL)
                CLEAN_AND_QUIT("main:getPath", "Failed to convert the input "
                        "to a string", SIFT3D_FALSE);

        // Convert the image to the internal format
        if (mx2imWithUnits(mxIm, mxUnits, &im))
                CLEAN_AND_QUIT("main:mx2im", "Failed to convert the input "
                        "image to the internal format", SIFT3D_TRUE);

        // Write the image
        switch (im_write(path, &im)) {
                case SIFT3D_SUCCESS:
                        break;
                case SIFT3D_UNSUPPORTED_FILE_TYPE:
                        CLEAN_AND_QUIT("main:unsupportedType", "Unsupported file "
                                "type", SIFT3D_TRUE);
                case SIFT3D_WRAPPER_NOT_COMPILED:
                        CLEAN_AND_QUIT("main:notCompiled", "Recompile SIFT3D "
                                "with support for this file type", SIFT3D_TRUE);
                default:
                        CLEAN_AND_QUIT("main:unexpected", "Unexpected error "
                                "writing the image", SIFT3D_TRUE);
        }

        // Clean up
        im_free(&im);

#undef CLEAN_AND_QUIT
}

