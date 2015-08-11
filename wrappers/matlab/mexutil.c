/* -----------------------------------------------------------------------------
 * mexutil.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex utility library for SIFT3D.
 * -----------------------------------------------------------------------------
 */

#include <stddef.h>
#include "mex.h"
#include "types.h"
#include "macros.h"
#include "imutil.h"
#include "sift.h"
#include "mexutil.h"

/* Keypoint struct information */
#define COORDS_NAME "coords"
#define SCALE_NAME "scale"
#define ORI_NAME "ori"
#define OCTAVE_NAME "octave"
#define LEVEL_NAME "level"
const char *fieldNames[] = {
        COORDS_NAME,
        SCALE_NAME,
        ORI_NAME,
        OCTAVE_NAME,
        LEVEL_NAME 
};
const mwSize kpNDims = 1;
const int kpNFields = sizeof(fieldNames) / sizeof(char *);

/* Entry point. */
SIFT3D sift3d;

/* Error message tag */
const char *tag = "sift3D";

/* Static helper functions */
static void init(void) __attribute__((constructor));
static void fini(void) __attribute__((destructor));

/* Library initialization */
static void init(void) {
        if (init_SIFT3D(&sift3d))
                err_msgu("main:initSift", "Failed to initialize SIFT3D");
}

/* Library cleanup */
static void fini(void) {
        cleanup_SIFT3D(&sift3d);
}

/* Print an error message. */
void err_msg(const char *name, const char *msg) {

        char id[1024];

        sprintf(id, "%s:%s", tag, name);

        mexErrMsgIdAndTxt(id, msg);
}

/* Print an unexpected error. */
void err_msgu(const char *name, const char *msg) {
        err_msg(name, msg);
        print_bug_msg();
}

/* Convert an mxArray to an Image. mx must have IM_NDIMS dimensions and be of
 * type single (float). */
int mx2im(const mxArray *const mx, Image *const im) {

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
mxArray *mat2mx(const Mat_rm *const mat) {

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

/* Convert from a Keypoint_store to an array of MATLAB keypoint structs. 
 * Returns the array of keypoints, or NULL on failure. */
mxArray *kp2mx(const Keypoint_store *const kp) {

        mxArray *mxKp;
        int i, coordsNum, scaleNum, oriNum, octaveNum, levelNum;

        const mwSize numKp = (mwSize) kp->slab.num;

	// Make an array of structs for the output 
	mxKp = mxCreateStructArray(kpNDims, &numKp, kpNFields, 
                                            fieldNames);
        if (mxKp == NULL)
                return NULL;

        // Get the field indices in the structs
        if ((coordsNum = mxGetFieldNumber(mxKp, COORDS_NAME)) < 0 ||
                (scaleNum = mxGetFieldNumber(mxKp, SCALE_NAME)) < 0 ||
                (oriNum = mxGetFieldNumber(mxKp, ORI_NAME)) < 0 ||
                (octaveNum = mxGetFieldNumber(mxKp, OCTAVE_NAME)) < 0 || 
                (levelNum = mxGetFieldNumber(mxKp, LEVEL_NAME)) < 0)
                return NULL;

        // Write the keypoints to the output
        for (i = 0; i < kp->slab.num; i++) {

                mxArray *mxCoords, *mxScale, *mxOri, *mxOctave, *mxLevel;
                double *coords;

                const Keypoint *const key = kp->buf + i;

                // Initialize the coordinate array
                if ((mxCoords = 
                        mxCreateDoubleMatrix(1, IM_NDIMS, mxREAL)) == NULL)
                        return NULL;

                // Copy the coordinates 
                coords = mxGetData(mxCoords); 
                coords[0] = key->xd;
                coords[1] = key->yd;
                coords[2] = key->zd;

                // Copy the transposed orientation
                if ((mxOri = mat2mx(&key->R)) == NULL)
                        return NULL;

                // Copy the scale 
                mxScale = mxCreateDoubleScalar(key->sd);

                // Copy the octave index
                mxOctave = mxCreateDoubleScalar((double) key->o); 

                // Copy the level index
                mxLevel = mxCreateDoubleScalar((double) key->s);
                
                // Set the struct fields
                mxSetFieldByNumber(mxKp, i, coordsNum, mxCoords);
                mxSetFieldByNumber(mxKp, i, scaleNum, mxScale);
                mxSetFieldByNumber(mxKp, i, oriNum, mxOri);
                mxSetFieldByNumber(mxKp, i, octaveNum, mxOctave);
                mxSetFieldByNumber(mxKp, i, levelNum, mxLevel);
        }

        return mxKp;
}

/* Convert an array of MATLAB keypoint structs to a Keypoint_store. */
int mx2kp(mxArray *mx, Keypoint_store *const kp) {

        //TODO

}

/* Wrapper for SIFT3D_detect_keypoints. */
int mex_SIFT3D_detect_keypoints(const Image *const im, 
        Keypoint_store *const kp) {
        return SIFT3D_detect_keypoints(&sift3d, im, kp);
}

/* Wrapper for SIFT3D_extract_descriptors. */
int mex_SIFT3D_extract_descriptors(const Image *const im, 
        const Keypoint_store *const kp, SIFT3D_Descriptor_store *const desc, 
        const int useGpyr) {
        return SIFT3D_extract_desciptors(&sift3d, im, &kp, &desc, useGpyr);
}
