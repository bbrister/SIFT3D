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

const char tag = "sift3D:main";

/* Keypoint struct information */
const char *coordsName = "coords";
const char *scaleName = "scale";
const char *oriName = "ori";
const char *fieldNames[] = {
        coordsName,
        scaleName,
        oriName 
};
const mwSize kpNDims = 1;
const int nfields = sizeof(fieldNames) / sizeof(char *);

/* Print an error message. */
static void err_msg(const char *name, const char *msg) {

        char id[1024];

        sprintf(id, "%s:%s", tag, name);

        mexErrMsgIdAndTxt(id, msg);
}

/* Print an unexpected error. */
static void err_msgu(const char *name, const char *msg) {
        err_msg(name, msg);
        print_bug_msg();
}

/* Convert a Mat_rm struct to an mxArray. Both must be initialized to the
 * correct size prior to calling this function. */
int mat2mxArray(const Mat_rm *const mat, mxArray *const mx) {

        double *mxData;
        int i, j;

        // TODO: get type of mx, convert mat to this type

        // Get the data
        if ((mxData = mxGetData(mx)) == NULL)
                return SIFT3D_FAILURE;

        // Copy a transposed version
        SIFT3D_MAT_RM_LOOP_START(mat, i, j)
                mxData[j + i * mat->num_cols] =
                        SIFT3D_MAT_RM_GET(mat, i, j, double);
        SIFT3D_MAT_RM_LOOP_END

        return SIFT3D_SUCCESS;
}

/* Entry point. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

        SIFT3D sift3d;
        mxArray *mxIm, *mxKp, *mxCoords, *mxScale, *mxOri;
        mwSize *imNDims;
        mwSize numKp;
        Image im;
        Keypoint_store kp;
        const char *errName, *errMsg;

/* Clean up and print an error */
#define CLEAN_AND_QUIT(name, msg) { \
                im_free(&im); \
                cleanup_Keypoint_store(&kp); \
                cleanup_SIFT3D(&sift3d); \
                err_msgu(name, msg); \
                if (mxCoords != NULL) \
                        mxDestroyArray(mxCoords); \
                if (mxScale != NULL) \
                        mxDestroyArray(mxScale); \
                if (mxOri != NULL) \
                        mxDestroyArray(mxOri); \
        }

	// Verify the number of inputs
        //TODO: Add another input for sift3d parameters
	if (nrhs > 1)
                err_msgu("numInputs", "This function takes 1 input.");

        // Verify the number of outputs
        if (nlhs != 1) 
                err_msgu("numOutputs", "This function takes 1 output.");

        // Verify that the image is of type float
	if (!mxIsSingle(mxIm) || mxIsComplex(mxIm)) 
                err_msgu("imType", "im must have type single");

	// Verify the number of image dimensions
	imNDims = mxGetNumberOfDimensions(mxIm);
	if (imNDims != 3) {
                err_msgu("imType", "im must have 3 dimensions");
        }

        // Assign the inputs
        mxIm = (mxArray *) prhs[0];

        // Initialize intermediates
        mxCoords = mxScale = mxOri = NULL;
        if (init_SIFT3D(&sift3d))
                err_msgu("initSift", "Failed to initialize SIFT3D");
        init_Keypoint_store(&kp);
        init_im(&im);
        if ((mxCoords = mxCreateDoubleMatrix(1, IM_NDIMS, mxREAL)) == NULL ||
                (mxScale = mxCreateDoubleMatrix(1, 1, mxScale)) == NULL ||
                (mxOri = mxCreateDoubleMatrix(IM_NDIMS, IM_NDIMS, mxScale)) 
                        == NULL)
                CLEAN_AND_QUIT("initTemp", "Failed to initialize temporaries");

        //TODO: parse the sift3d parameters 

        // Detect keypoints
	if (SIFT3D_detect_keypoints(&sift3d, &im, &kp))
                CLEAN_AND_QUIT("detect", "Failed to detect keypoints");

	// Make an array of structs for the output 
	mxKp = plhs[0] = mxCreateStructArray(kpNDims, &numKp, nFields, 
                                            fieldNames);
        if (mxKp == NULL)
                CLEAN_AND_QUIT("createOutput", "Failed to create outputs");


        // Write the keypoints to the output
        for (int i = 0; i < kp.slab.num; i++) {

                const Keypoint *const key = kp->buf + i;

                // Copy the coordinates 
                coords[0] = key->xd;
                coords[1] = key->yd;
                coords[2] = key->zd;

                // Copy the scale
                scale[0] = key->sd; 

                // Copy the transposed orientation
                if (mat2mxArray(&kp->ori, mxOri))
                        CLEAN_AND_QUIT("convertOri", "Failed to convert "
                                "orientation matrix");

                // Set the struct fields
                mxSetField(mxKp, i, coordsName, mxCoords);
                mxSetField(mxKp, i, scaleName, mxScale);
                mxSetField(mxKp, i, oriName, mxOri);
        }

        // Clean up
        im_free(&im);
        cleanup_Keypoint_store(&kp);
        cleanup_SIFT3D(&sift3d);
        if (mxCoords != NULL)
                mxDestroyArray(mxCoords);
        if (mxScale != NULL)
                mxDestroyArray(mxScale);
        if (mxOri != NULL)
                mxDestroyArray(mxOri);

#undef CLEAN_AND_QUIT
}
