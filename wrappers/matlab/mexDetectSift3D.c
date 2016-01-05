/* -----------------------------------------------------------------------------
 * mexDetectSift3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
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

        const mxArray *mxIm, *mxUnits;
        Image im;
        Keypoint_store kp;
        const char *errName, *errMsg;

/* Clean up and print an error */
#define CLEAN_AND_QUIT(name, msg, expected) { \
                im_free(&im); \
                cleanup_Keypoint_store(&kp); \
                if (expected) { \
                        err_msg(name, msg); \
                } else { \
                        err_msgu(name, msg); \
                } \
        }

	// Verify the number of inputs
	if (nrhs != 2)
                err_msg("main:numInputs", "This function takes 2 inputs.");

        // Verify the number of outputs
        if (nlhs > 1) 
                err_msg("main:numOutputs", "This function takes 1 output.");

        // Assign the inputs
        mxIm = prhs[0];
        mxUnits = prhs[1];

        // Initialize intermediates
        init_Keypoint_store(&kp);
        init_im(&im);

        // Copy the transposed image
        if (mx2imWithUnits(mxIm, mxUnits, &im))
                CLEAN_AND_QUIT("main:copyIm", "Failed to convert the input "
                        "image", SIFT3D_TRUE);

        // Detect keypoints
	if (mex_SIFT3D_detect_keypoints(&im, &kp))
                CLEAN_AND_QUIT("main:detect", "Failed to detect keypoints", 
                        SIFT3D_FALSE);

        // Convert the output to a MATLAB array of structs
        if ((plhs[0] = kp2mx(&kp)) == NULL)
                CLEAN_AND_QUIT("main:convertToStructs", "Failed to convert "
                        "keypoints to structs", SIFT3D_FALSE);

        // Clean up
        im_free(&im);
        cleanup_Keypoint_store(&kp);

#undef CLEAN_AND_QUIT
}

