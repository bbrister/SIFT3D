/* -----------------------------------------------------------------------------
 * mexImWrite3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
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
#define CLEAN_AND_QUIT(name, msg) { \
                im_free(&im); \
                err_msgu(name, msg); \
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
                        "to a string");

        // Convert the image to the internal format
        if (mx2im(mxIm, &im))
                CLEAN_AND_QUIT("main:mx2im", "Failed to convert the input "
                        "image array to the internal image format");

        // Convert the units to the internal format
        if (mx2units(mxUnits, &im))
                CLEAN_AND_QUIT("main:mx2im", "Failed to convert the input "
                        "units array to the internal format");

        // Write the image
        switch (im_write(path, &im)) {
                case SIFT3D_SUCCESS:
                        break;
                case SIFT3D_UNSUPPORTED_FILE_TYPE:
                        im_free(&im);
                        err_msg("main:unsupportedType", "Unsupported file "
                                "type");
                        break;
                default:
                        CLEAN_AND_QUIT("main:unexpected", "Unexpected error "
                                "reading the image");
        }

        // Clean up
        im_free(&im);

#undef CLEAN_AND_QUIT
}

