/* -----------------------------------------------------------------------------
 * mexDetectSift3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex interface to detect SIFT3D keypoints.
 * -----------------------------------------------------------------------------
 */

#include "mex.h"
#include "sift.h"
#include "imutil.h"
#include "macros.h"
#include "mexutil.h"

/* Keypoint struct information */
#define COORDS_NAME "coords"
#define SCALE_NAME "scale"
#define ORI_NAME "ori"
const char *fieldNames[] = {
        COORDS_NAME,
        SCALE_NAME,
        ORI_NAME 
};
const mwSize kpNDims = 1;
const int nFields = sizeof(fieldNames) / sizeof(char *);

/* Entry point. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

        const mxArray *mxIm;
        mxArray *mxKp;
        mwSize imNDims, numKp;
        Image im;
        Keypoint_store kp;
        const char *errName, *errMsg;
        int i, coordsNum, scaleNum, oriNum;

/* Clean up and print an error */
#define CLEAN_AND_QUIT(name, msg) { \
                im_free(&im); \
                cleanup_Keypoint_store(&kp); \
                err_msgu(name, msg); \
        }

	// Verify the number of inputs
	if (nrhs > 1)
                err_msgu("main:numInputs", "This function takes 1 input.");

        // Verify the number of outputs
        if (nlhs != 1) 
                err_msgu("main:numOutputs", "This function takes 1 output.");

        // Assign the inputs
        mxIm = prhs[0];

        // Verify that the image is of type float
	if (!mxIsSingle(mxIm) || mxIsComplex(mxIm)) 
                err_msgu("main:imType", "im must have type single");

	// Verify the number of image dimensions
	imNDims = mxGetNumberOfDimensions(mxIm);
	if (imNDims != IM_NDIMS) {

                char msg[1024];

                sprintf(msg, "im must have %d dimensions", IM_NDIMS);
                err_msgu("main:imType", msg);
        }

        // Initialize intermediates
        init_Keypoint_store(&kp);
        init_im(&im);

        // Copy the transposed image
        if (mx2im(mxIm, &im))
                CLEAN_AND_QUIT("main:copyIm", "Failed to copy the input image");

        // Detect keypoints
	if (mex_SIFT3D_detect_keypoints(&im, &kp))
                CLEAN_AND_QUIT("main:detect", "Failed to detect keypoints");

	// Make an array of structs for the output 
        numKp = (mwSize) kp.slab.num;
	mxKp = plhs[0] = mxCreateStructArray(kpNDims, &numKp, nFields, 
                                            fieldNames);
        if (mxKp == NULL)
                CLEAN_AND_QUIT("main:createOutput", "Failed to create outputs");

        // Get the field indices in the structs
        if ((coordsNum = mxGetFieldNumber(mxKp, COORDS_NAME)) < 0 ||
                (scaleNum = mxGetFieldNumber(mxKp, SCALE_NAME)) < 0 ||
                (oriNum = mxGetFieldNumber(mxKp, ORI_NAME)) < 0)
                CLEAN_AND_QUIT("main:getFields", "Failed to get field indices");

        // Write the keypoints to the output
        for (i = 0; i < kp.slab.num; i++) {

                mxArray *mxCoords, *mxScale, *mxOri;
                double *coords, *scale;

                const Keypoint *const key = kp.buf + i;

                // Initialize arrays for the data
                if ((mxCoords = 
                        mxCreateDoubleMatrix(1, IM_NDIMS, mxREAL)) == NULL ||
                        (mxScale =
                                mxCreateDoubleMatrix(1, 1, mxREAL)) == NULL)
                        CLEAN_AND_QUIT("main:initArrays", 
                                "Failed to initialize field arrays");

                // Get pointers to the internal data of the field matrices
                coords = mxGetData(mxCoords); 
                scale = mxGetData(mxScale); 

                // Copy the coordinates 
                coords[0] = key->xd;
                coords[1] = key->yd;
                coords[2] = key->zd;

                // Copy the scale
                scale[0] = key->sd; 
        
                // Copy the transposed orientation
                if ((mxOri = mat2mx(&key->R)) == NULL)
                        CLEAN_AND_QUIT("main:convertOri", "Failed to convert "
                                "orientation matrix");

                // Set the struct fields
                mxSetFieldByNumber(mxKp, i, coordsNum, mxCoords);
                mxSetFieldByNumber(mxKp, i, scaleNum, mxScale);
                mxSetFieldByNumber(mxKp, i, oriNum, mxOri);
        }

        // Clean up
        im_free(&im);
        cleanup_Keypoint_store(&kp);

#undef CLEAN_AND_QUIT
}

