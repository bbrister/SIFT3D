/* -----------------------------------------------------------------------------
 * mexExtractSift3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex interface to extract SIFT3D descriptors.
 * -----------------------------------------------------------------------------
 */

#include "mex.h"
#include "matrix.h"
#include "imutil.h"
#include "macros.h"
#include "mexutil.h"

/* Constants */
#define DESC_NDIMS 2

/* Entry point. 
 * Output format:
 * [x y z el0 el1 ... el767] (see sift.c:SIFT3D_Descriptor_store_to_Mat_rm)
 *
 * Note that the matlab function does some postprocessing to present the data
 * in a different output format.
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

        const mxArray *mxKp, *mxIm;
        const char *errMsg;
        const void *imArg;
        Image im;
        Keypoint_store kp;
        SIFT3D_Descriptor_store desc;
        int i, useGpyr;

/* Clean up and print an error */
#define CLEAN_AND_QUIT(name, msg) { \
                im_free(&im); \
                cleanup_Keypoint_store(&kp); \
                cleanup_SIFT3D_Descriptor_store(&desc); \
                err_msgu(name, msg); \
        }

	// Verify the number of inputs
	if (nrhs != 2)
                err_msgu("main:numInputs", "This function takes 2 inputs.");

        // Verify the number of outputs
        if (nlhs > 1) 
                err_msgu("main:numOutputs", "This function takes 1 output.");

        // Assign inputs
        mxKp = prhs[0];
        mxIm = prhs[1];

        // Initialize intermediates
        init_im(&im); 
        init_Keypoint_store(&kp);
        init_SIFT3D_Descriptor_store(&desc);

        // Convert the keypoints
        if (mx2kp(mxKp, &kp))
                CLEAN_AND_QUIT("main:convertKp", "Failed to convert keypoints");

        // Convert the image, if any
        if (!mxIsEmpty(mxIm)) {
                if (mx2im(mxIm, &im))
                        CLEAN_AND_QUIT("main:convertIm", 
                                        "Failed to convert image");
                imArg = &im;
                useGpyr = SIFT3D_FALSE;
        } else {
                if ((imArg = mexGetGpyr()) == NULL)
                        CLEAN_AND_QUIT("main:getGpyr",
                                "Failed to get the gaussian pyramid");
                useGpyr = SIFT3D_TRUE; 
        }

        // Extract descriptors
        if (mex_SIFT3D_extract_descriptors(imArg, &kp, &desc, useGpyr))
                err_msgu("main:extract", "Failed to extract descriptors");

        // Convert the descriptors to an output matrix
        if ((plhs[0] = desc2mx(&desc)) == NULL)
                CLEAN_AND_QUIT("main:createOutput", 
                        "Failed to convert descriptors");

        // Clean up
        im_free(&im);
        cleanup_Keypoint_store(&kp);
        cleanup_SIFT3D_Descriptor_store(&desc);

#undef CLEAN_AND_QUIT
}
