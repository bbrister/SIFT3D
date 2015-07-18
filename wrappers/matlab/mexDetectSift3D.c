/* -----------------------------------------------------------------------------
 * mexDetectSift3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex interface to detect SIFT3D keypoints.
 * -----------------------------------------------------------------------------
 */

#include <stddef.h>
#include "mex.h"
#include "sift.h"
#include "imutil.h"
#include "macros.h"

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

/* Error message tag */
const char *tag = "sift3D";

/* Helper function declarations. */
static void err_msg(const char *name, const char *msg);
static void err_msgu(const char *name, const char *msg);
static mxArray *mat2mxArray(const Mat_rm *const mat, mxArray *const mx);

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

/* Convert an mxArray to an Image. mx must have IM_NDIMS dimensions and be of
 * type single (float). */
static int mx2im(const mxArray *const mx, Image *const im) {

        const mwSize *mxDims;
        float *mxData;
        int i, x, y, z;

        // Verify inputs
	if (mxGetNumberOfDimensions(mx) != IM_NDIMS ||
                !mxIsSingle(mx))
                return SIFT3D_FAILURE;

        // Copy the dimensions
        mxDims = mxGetDimensions(mx);
        for (i = 0; i < IM_NDIMS; i++) {
                im->dims[i] = (int) mxDims[i];
        }

        // Resize the output                
        im->nc = 1;
        im_default_stride(im);
        if (im_resize(im))
                return SIFT3D_FAILURE;

        // Get the data
        if ((mxData = mxGetData(mx)) == NULL)
                return SIFT3D_FAILURE;

        // Transpose and copy the data
        SIFT3D_IM_LOOP_START(im, x, y, z)

                mwIndex idx;

                const mwIndex subs[] = {x, y, z};

                idx = mxCalcSingleSubscript(mx, *mxDims, subs);
                SIFT3D_IM_GET_VOX(im, x, y, z, 0) = mxData[idx];

        SIFT3D_IM_LOOP_END

        return SIFT3D_SUCCESS;
}

/* Convert a Mat_rm struct to an mxArray. Returns the array, or NULL on
 * failure. The returned array has type double, regardless of the input array
 * type. */
static mxArray *mat2mx(const Mat_rm *const mat) {

        mxArray *mx;
        double *mxData;
        int i, j;

        const int rows = mat->num_rows;
        const int cols = mat->num_cols;

        // Create an array
        if ((mx = mxCreateDoubleMatrix(rows, cols, mxREAL)) == NULL)
                return NULL;

        // Get the data
        if ((mxData = mxGetData(mx)) == NULL)
                goto mat2mx_quit;

#define SIFT3D_MAT_TYPE_MACRO(mat, err_label, MACRO, ...) \
        switch ((mat)->type) { \
                case DOUBLE: \
                        MACRO(double, ## __VA_ARGS__) \
                        break; \
                case FLOAT: \
                        MACRO(float, ## __VA_ARGS__) \
                        break; \
                case INT: \
                        MACRO(int, ## __VA_ARGS__) \
                        break; \
                default: \
                        fprintf(stderr, "imutil: unknown matrix type"); \
                        goto err_label; \
        } \

#define TRANSPOSE_AND_COPY(type) \
        SIFT3D_MAT_RM_LOOP_START(mat, i, j) \
        \
                mwIndex idx; \
        \
                const mwIndex subs[] = {i, j}; \
        \
                idx = mxCalcSingleSubscript(mx, 2, subs); \
                mxData[idx] = (double) SIFT3D_MAT_RM_GET(mat, i, j, type); \
        \
        SIFT3D_MAT_RM_LOOP_END

        // Transpose and copy the data 
        SIFT3D_MAT_TYPE_MACRO(mat, mat2mx_quit, TRANSPOSE_AND_COPY);

#undef TRANSPOSE_AND_COPY

        return mx;

mat2mx_quit:
        mxDestroyArray(mx);
        return NULL;
}

/* Entry point. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

        SIFT3D sift3d;
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
                cleanup_SIFT3D(&sift3d); \
                err_msgu(name, msg); \
        }

	// Verify the number of inputs
        //TODO: Add another input for sift3d parameters
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
        if (init_SIFT3D(&sift3d))
                err_msgu("main:initSift", "Failed to initialize SIFT3D");
        init_Keypoint_store(&kp);
        init_im(&im);

        //TODO: parse the sift3d parameters 

        // Copy the transposed image
        if (mx2im(mxIm, &im))
                CLEAN_AND_QUIT("main:copyIm", "Failed to copy the input image");

        // Detect keypoints
	if (SIFT3D_detect_keypoints(&sift3d, &im, &kp))
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
        cleanup_SIFT3D(&sift3d);

#undef CLEAN_AND_QUIT
}

