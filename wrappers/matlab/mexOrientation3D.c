/* -----------------------------------------------------------------------------
 * mexOrientation3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex interface to assign 3D orientation.
 * -----------------------------------------------------------------------------
 */

#include "imutil.h"
#include "sift.h"
#include "immacros.h"
#include "mexutil.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

        const mxArray *mxKp, *mxIm, *mxUnits;
        Image im;
        Keypoint_store kp;
        const char *errName, *errMsg;
        double *conf;

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
	if (nrhs != 3)
                err_msg("main:numInputs", "This function takes 3 inputs.");

        // Verify the number of outputs
        if (nlhs > 2) 
                err_msg("main:numOutputs", "This function takes 2 outputs.");

        // Assign the inputs
        mxKp = prhs[0];
        mxIm = prhs[1];
        mxUnits = prhs[2];

        // Initialize intermediates
        init_im(&im);
        init_Keypoint_store(&kp);

        // Convert the keypoints 
        if (mx2kp(mxKp, &kp))
                CLEAN_AND_QUIT("main:convertKey", "failed to convert the "
                        "input keypoints", SIFT3D_TRUE);

        // Convert the image
        if (mx2imWithUnits(mxIm, mxUnits, &im))
                CLEAN_AND_QUIT("main:copyIm", "Failed to convert the input "
                        "image", SIFT3D_TRUE);

        // Assign the orientations
        conf = NULL;
        if (mex_SIFT3D_assign_orientations(&im, &kp, &conf))
                CLEAN_AND_QUIT("main:assignOrientations", "Failed to assign "
                        "the orientations", SIFT3D_TRUE);

        // Convert the output keypoints
        if ((plhs[0] = kp2mx(&kp)) == NULL)
                CLEAN_AND_QUIT("main:convertR", "Failed to convert R", 
                        SIFT3D_FALSE);

        // Convert the output confidence to a MATLAB array
        if ((plhs[1] = array2mx(conf, kp.slab.num)) == NULL)
                CLEAN_AND_QUIT("main:convertConf", "Failed to convert conf",
                        SIFT3D_FALSE);

        // Clean up
        im_free(&im);
        cleanup_Keypoint_store(&kp);

#undef CLEAN_AND_QUIT
}

