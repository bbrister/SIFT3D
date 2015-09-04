/* -----------------------------------------------------------------------------
 * mexDetectSift3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex interface to detect SIFT3D keypoints.
 * -----------------------------------------------------------------------------
 */

#include "imutil.h"
#include "sift.h"
#include "macros.h"
#include "mexutil.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

        const mxArray *mxIm;
        mwSize imNDims;
        Image im;
        Keypoint_store kp;
        const char *errName, *errMsg;

/* Clean up and print an error */
#define CLEAN_AND_QUIT(name, msg) { \
                im_free(&im); \
                cleanup_Keypoint_store(&kp); \
                err_msgu(name, msg); \
        }

	// Verify the number of inputs
	if (nrhs != 1)
                err_msgu("main:numInputs", "This function takes 1 input.");

        // Verify the number of outputs
        if (nlhs > 1) 
                err_msgu("main:numOutputs", "This function takes 1 output.");

        // Assign the inputs
        mxIm = prhs[0];

        // Initialize intermediates
        init_Keypoint_store(&kp);
        init_im(&im);

        // Copy the transposed image
        if (mx2im(mxIm, &im))
                CLEAN_AND_QUIT("main:copyIm", "Failed to copy the input image");

        // Detect keypoints
	if (mex_SIFT3D_detect_keypoints(&im, &kp))
                CLEAN_AND_QUIT("main:detect", "Failed to detect keypoints");

        // Convert the output to a MATLAB array of structs
        if ((plhs[0] = kp2mx(&kp)) == NULL)
                CLEAN_AND_QUIT("main:convertToStructs", "Failed to convert "
                        "keypoints to structs");

        // Clean up
        im_free(&im);
        cleanup_Keypoint_store(&kp);

#undef CLEAN_AND_QUIT
}

