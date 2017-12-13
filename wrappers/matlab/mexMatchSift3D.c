/* -----------------------------------------------------------------------------
 * mexMatchSift3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2017 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex interface to register images from SIFT3D keypoints and descriptors.
 * -----------------------------------------------------------------------------
 */

#include <stdlib.h>
#include "imutil.h"
#include "immacros.h"
#include "sift.h"
#include "reg.h"
#include "mexutil.h"
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

        const mxArray *mxDesc1, *mxDesc2, *mxNnThresh;
        const mwSize *mxDims;
        mxArray *mxMatches;
        double *mxMatchesData;
        int *matches;
        SIFT3D_Descriptor_store desc1, desc2;
        double nn_thresh;
        int i;

/* Clean up and print an error */
#define CLEAN_AND_QUIT(name, msg, expected) { \
                cleanup_SIFT3D_Descriptor_store(&desc1); \
                cleanup_SIFT3D_Descriptor_store(&desc2); \
                if (matches != NULL) { \
                        free(matches);  \
                } \
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
        if (nlhs > 1)
                err_msg("main:numOutputs", "This function takes 1 output.");

        // Assign the inputs
        mxDesc1 = prhs[0];
        mxDesc2 = prhs[1];
        mxNnThresh = prhs[2];

        // Get the matching threshold, if one was provided
        nn_thresh = mxIsEmpty(mxNnThresh) ? SIFT3D_nn_thresh_default :
                mxGetScalar(mxNnThresh);

        // Initialize intermediates
        matches = NULL;
	init_SIFT3D_Descriptor_store(&desc1);
	init_SIFT3D_Descriptor_store(&desc2);

        // Convert the inputs to descriptor stores
        if (mx2desc(mxDesc1, &desc1))
                CLEAN_AND_QUIT("main:convert1", "Failed to convert desc1",
                        SIFT3D_TRUE);
        if (mx2desc(mxDesc2, &desc2))
                CLEAN_AND_QUIT("main:convert2", "Failed to convert desc2",
                        SIFT3D_TRUE);

        // Match descriptors
	if (SIFT3D_nn_match(&desc1, &desc2, nn_thresh, &matches))
                CLEAN_AND_QUIT("main:match", "Failed to match descriptors",
                        SIFT3D_FALSE);

        // Create a matrix for the output
        if ((mxMatches = mxCreateDoubleMatrix(desc1.num, 1, mxREAL)) == NULL)
                CLEAN_AND_QUIT("main:createOutput", "Failed to create output",
                        SIFT3D_TRUE);
        plhs[0] = mxMatches;

        // Get the data
        if ((mxMatchesData = mxGetData(mxMatches)) == NULL)
                CLEAN_AND_QUIT("main:convertOut", "Failed to get output data",
                        SIFT3D_FALSE);

        // Copy the output
        for (i = 0; i < desc1.num; i++) {
                mxMatchesData[i] = (double) matches[i];
        }

        // Clean up
        cleanup_SIFT3D_Descriptor_store(&desc1);
        cleanup_SIFT3D_Descriptor_store(&desc2);
        free(matches);
}
