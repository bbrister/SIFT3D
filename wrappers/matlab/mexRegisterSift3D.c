/* -----------------------------------------------------------------------------
 * mexRegisterSift3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex interface to register images from SIFT3D keypoints and descriptors.
 * -----------------------------------------------------------------------------
 */

#include "imutil.h"
#include "reg.h"
#include "macros.h"
#include "mexutil.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

        const mxArray *mxSrc, *mxRef, *mxSrcUnits, *mxRefUnits, *mxThresh;
        Image src, ref;
        Affine aff;
        Mat_rm match_src, match_ref;
        double nn_thresh_old;
        int ret;

/* Clean up and print an error */
#define CLEAN_AND_QUIT(name, msg, expected) { \
                im_free(&src); \
                im_free(&ref); \
                cleanup_tform(&aff); \
                cleanup_Mat_rm(&match_src); \
                cleanup_Mat_rm(&match_ref); \
                if (expected) { \
                        err_msg(name, msg); \
                } else { \
                        err_msgu(name, msg); \
                } \
        }

	// Verify the number of inputs
	if (nrhs != 5)
                err_msg("main:numInputs", "This function takes 5 inputs.");

        // Verify the number of outputs
        if (nlhs > 3)
                err_msg("main:numOutputs", "This function takes 3 outputs.");

        // Assign the inputs
        mxSrc = prhs[0];
        mxRef = prhs[1];
        mxSrcUnits = prhs[2];
        mxRefUnits = prhs[3];
        mxThresh = prhs[4];

        // Initialize intermediates
        if (init_Affine(&aff, IM_NDIMS) ||
                init_Mat_rm(&match_src, 0, 0, DOUBLE, SIFT3D_FALSE) ||
                init_Mat_rm(&match_ref, 0, 0, DOUBLE, SIFT3D_FALSE))
                err_msgu("main:init", "Failed to initialize intermediates");
        init_im(&src);
        init_im(&ref);

        // Convert the inputs to images
        if (mx2imWithUnits(mxSrc, mxSrcUnits, &src))
                CLEAN_AND_QUIT("main:convertSrc", "Failed to convert the "
                        "source image.", SIFT3D_TRUE);
        if (mx2imWithUnits(mxRef, mxRefUnits, &ref))
                CLEAN_AND_QUIT("main:convertSrc", "Failed to convert the "
                        "reference image.", SIFT3D_TRUE);
        
        // Set the threshold if it is not empty
        if (!mxIsEmpty(mxThresh)) {

                // Verify thresh
                if (!mxIsScalar(mxThresh)) {
                        CLEAN_AND_QUIT("main:threshScalar", 
                                "thresh must be scalar.", SIFT3D_TRUE);
                }
        
                // Record the old value
                nn_thresh_old = mex_get_nn_thresh_Reg_SIFT3D();
        
                // Set the new one
                if (mex_set_nn_thresh_Reg_SIFT3D(mxGetScalar(mxThresh)))
                        CLEAN_AND_QUIT("main:thresh", "Invalid thresh.", 
                                SIFT3D_TRUE);
        }

        // Register the images
        ret = mex_register_SIFT3D_resample(&src, &ref, LINEAR, &aff);

        // Optionally reset the old threshold
        if (!mxIsEmpty(mxThresh) && 
                mex_set_nn_thresh_Reg_SIFT3D(nn_thresh_old))
                CLEAN_AND_QUIT("main:resetThresh", "Failed to reset thresh.", 
                        SIFT3D_FALSE);

        // Handle registration errors
        if (ret) {
                CLEAN_AND_QUIT("main:reg", "Failed to register the images. "
                        "Possibly no good transformation was found.", 
                        SIFT3D_TRUE);
        }

        // Retrieve the matches
        if (mex_get_matches_Reg_SIFT3D(&match_src, &match_ref))
                CLEAN_AND_QUIT("main:getMatches", "Failed to retrieve the "
                        "matches.", SIFT3D_FALSE);

        // Convert the outputs to mxArrays
        if ((plhs[0] = mat2mx(&aff.A)) == NULL ||
                (plhs[1] = mat2mx(&match_src)) == NULL ||
                (plhs[2] = mat2mx(&match_ref)) == NULL)
                CLEAN_AND_QUIT("main:convertA", "Failed to convert "
                        "outputs.", SIFT3D_FALSE);

        // Clean up
        im_free(&src);
        im_free(&ref); 
        cleanup_tform(&aff); 
        cleanup_Mat_rm(&match_src); 
        cleanup_Mat_rm(&match_ref); 

#undef CLEAN_AND_QUIT        
}

