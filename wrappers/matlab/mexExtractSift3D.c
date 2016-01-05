/* -----------------------------------------------------------------------------
 * mexExtractSift3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex interface to extract SIFT3D descriptors.
 * -----------------------------------------------------------------------------
 */

#include "imutil.h"
#include "sift.h"
#include "mexutil.h"
#include "mex.h"

/* Entry point. 
 * Output format:
 * [x y z el0 el1 ... el767] (see sift.c:SIFT3D_Descriptor_store_to_Mat_rm)
 *
 * Note that the matlab function does some postprocessing to present the data
 * in a different output format.
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

        const mxArray *mxKp, *mxIm, *mxUnits;
        const char *errMsg;
        Image im;
        Keypoint_store kp;
        SIFT3D_Descriptor_store desc;
        int i;

/* Clean up and print an error */
#define CLEAN_AND_QUIT(name, msg, expected) { \
                im_free(&im); \
                cleanup_Keypoint_store(&kp); \
                cleanup_SIFT3D_Descriptor_store(&desc); \
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
        if (nlhs > 1) 
                err_msgu("main:numOutputs", "This function takes 1 output.");

        // Assign inputs
        mxKp = prhs[0];
        mxIm = prhs[1];
        mxUnits = prhs[2];

        // Initialize intermediates
        init_im(&im); 
        init_Keypoint_store(&kp);
        init_SIFT3D_Descriptor_store(&desc);

        // Convert the keypoints
        if (mx2kp(mxKp, &kp))
                CLEAN_AND_QUIT("main:convertKp", "Failed to convert keypoints",
                        SIFT3D_TRUE);

        // Process the image and extract descriptors
        if (!mxIsEmpty(mxIm)) {

                // Convert the input to an Image struct
                if (mx2imWithUnits(mxIm, mxUnits, &im))
                        CLEAN_AND_QUIT("main:convertIm", 
                                        "Failed to convert image", 
                                        SIFT3D_TRUE);

                // Extract raw descriptors
                if (mex_SIFT3D_extract_raw_descriptors(&im, &kp, &desc))
                        CLEAN_AND_QUIT("main:extractRaw", 
                                "Failed to extract raw descriptors",
                                SIFT3D_TRUE);
        } else {

                // Attempt to retrieve the Gaussian pyramid
                if (!mexHaveGpyr)
                        CLEAN_AND_QUIT("main:haveGpyr",
                                "Failed to get the Gaussian pyramid. Must "
                                "call detectSift3D before this function can "
                                "be called without the im argument.", 
                                SIFT3D_TRUE);

                // Extract descriptors from the pyramid
                if (mex_SIFT3D_extract_descriptors(&kp, &desc))
                        CLEAN_AND_QUIT("main:extractPyramid", 
                                "Failed to extract pyramid descriptors",
                                SIFT3D_FALSE);
        }


        // Convert the descriptors to an output matrix
        if ((plhs[0] = desc2mx(&desc)) == NULL)
                CLEAN_AND_QUIT("main:createOutput", 
                        "Failed to convert descriptors", SIFT3D_FALSE);

        // Clean up
        im_free(&im);
        cleanup_Keypoint_store(&kp);
        cleanup_SIFT3D_Descriptor_store(&desc);

#undef CLEAN_AND_QUIT
}
