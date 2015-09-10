/* -----------------------------------------------------------------------------
 * mexutil.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex utility library for SIFT3D.
 * -----------------------------------------------------------------------------
 */

/* Standard headers */
#include <stddef.h>
#include <math.h>

/* SIFT3D headers */
#include "types.h"
#include "macros.h"
#include "imutil.h"
#include "sift.h"
#include "mexutil.h"

/* Matlab headers */
#include "mex.h"
#include "matrix.h"

/* Execute the macro MACRO, with the first argument set to the type of mat. If
 * there is an error, goto err_label. */
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

/* Returns SIFT3D_TRUE if the input is a real-valued double-precision floating
 * point array, aka "double" in C. */
int isDouble(const mxArray *const mx) {
        return mxIsDouble(mx) && !mxIsComplex(mx);
}

/* Convert an mxArray to an Image. mx must have IM_NDIMS dimensions and be of
 * type single (float). */
int mx2im(const mxArray *const mx, Image *const im) {

        const mwSize *mxDims;
        float *mxData;
        int i, x, y, z;

        // Verify inputs
	if (mxGetNumberOfDimensions(mx) != IM_NDIMS ||
                !mxIsSingle(mx) || mxIsComplex(mx))
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

                idx = mxCalcSingleSubscript(mx, IM_NDIMS, subs);
                SIFT3D_IM_GET_VOX(im, x, y, z, 0) = mxData[idx];

        SIFT3D_IM_LOOP_END

        return SIFT3D_SUCCESS;
}

/* Convert a Mat_rm struct to an mxArray. Returns the array, or NULL on
 * failure. The returned array has type double, regardless of the input array
 * type. 
 *
 * Returns the array, or NULL on failure. */
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

/* Convert an mxArray to a Mat_rm. mx must have type double. The type of mat
 * is preserved.
 * 
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */ 
int mx2mat(const mxArray *const mx, Mat_rm *const mat) {

        const mwSize *mxDims;
        double *data;
        int i, j;

        // Verify inputs
        if (!isDouble(mx) || mxGetNumberOfDimensions(mx) != 2)
                return SIFT3D_FAILURE;

        // Resize the output
        mxDims = mxGetDimensions(mx);
        mat->num_rows = (int) mxDims[0];
        mat->num_cols = (int) mxDims[1];
        if (resize_Mat_rm(mat))
                return SIFT3D_FAILURE;

        // Get the data
        if ((data = mxGetData(mx)) == NULL)
                return SIFT3D_FAILURE;

#define COPY_DATA(type) \
        SIFT3D_MAT_RM_LOOP_START(mat, i, j) \
\
                mwIndex idx; \
\
                const mwIndex subs[] = {i, j}; \
\
                idx = mxCalcSingleSubscript(mx, 2, subs); \
                SIFT3D_MAT_RM_GET(mat, i, j, type) = (type) data[idx]; \
\
        SIFT3D_MAT_RM_LOOP_END

        // Copy the transposed data
        SIFT3D_MAT_TYPE_MACRO(mat, mx2mat_quit, COPY_DATA);

#undef COPY_DATA

        return SIFT3D_SUCCESS;

mx2mat_quit:
        return SIFT3D_FAILURE;
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

/* Convert an array of MATLAB keypoint structs to a Keypoint_store. mx must be
 * a vector, and each struct element must have type double. */
int mx2kp(const mxArray *const mx, Keypoint_store *const kp) {

        const mwSize *mxDims, *coordsDims, *oriDims;
        mxArray *mxCoords, *mxScale, *mxOri, *mxOctave, *mxLevel;
        mwSize numKp, kpRows, kpCols;
        int i, coordsNum, scaleNum, oriNum, octaveNum, levelNum;

        // Verify the number of dimensions
        if (mxGetNumberOfDimensions(mx) != 2)
                return SIFT3D_FAILURE;

        // Parse the input dimensions
        mxDims = mxGetDimensions(mx);
        kpRows = mxDims[0];
        kpCols = mxDims[1];
        if (kpRows == 1)
                numKp = kpCols; 
        else if (kpCols == 1)
                numKp = kpRows;
        else
                return SIFT3D_FAILURE; 

        // Verify the struct size 
        if (!mxIsStruct(mx) || mxGetNumberOfFields(mx) != kpNFields)
                return SIFT3D_FAILURE;

        // Get the field indices
        if ((coordsNum = mxGetFieldNumber(mx, COORDS_NAME)) < 0 ||
                (scaleNum = mxGetFieldNumber(mx, SCALE_NAME)) < 0 ||
                (oriNum = mxGetFieldNumber(mx, ORI_NAME)) < 0 ||
                (octaveNum = mxGetFieldNumber(mx, OCTAVE_NAME)) < 0 || 
                (levelNum = mxGetFieldNumber(mx, LEVEL_NAME)) < 0)
                return SIFT3D_FAILURE;

        // Get the fields of the first keypoint
        if ((mxCoords = mxGetFieldByNumber(mx, 0, coordsNum)) == NULL ||
                (mxScale = mxGetFieldByNumber(mx, 0, scaleNum)) == NULL ||
                (mxOri = mxGetFieldByNumber(mx, 0, oriNum)) == NULL ||
                (mxOctave = mxGetFieldByNumber(mx, 0, octaveNum)) == NULL ||
                (mxLevel = mxGetFieldByNumber(mx, 0, levelNum)) == NULL)
                return SIFT3D_FAILURE;

        // Verify the type
        if (!isDouble(mxCoords) ||
                !isDouble(mxScale) ||
                !isDouble(mxOri) ||
                !isDouble(mxOctave) ||
                !isDouble(mxLevel))
                return SIFT3D_FAILURE;

        // Verify the number of dimensions
        if (mxGetNumberOfDimensions(mxCoords) != 2 ||
                mxGetNumberOfDimensions(mxScale) != 2 ||
                mxGetNumberOfDimensions(mxOri) != 2 ||
                mxGetNumberOfDimensions(mxOctave) != 2 ||
                mxGetNumberOfDimensions(mxLevel) != 2)
                return SIFT3D_FAILURE;

        // Verify the scalars
        if (!mxIsScalar(mxScale) ||
                !mxIsScalar(mxOctave) ||
                !mxIsScalar(mxLevel))
                return SIFT3D_FAILURE;

        // Verify the coordinate vector dimensions 
        coordsDims = mxGetDimensions(mxCoords);
        if ((coordsDims[0] != 1 && coordsDims[1] != 1) ||
                coordsDims[0] * coordsDims[1] != IM_NDIMS)
                return SIFT3D_FAILURE;


        // Verify the orientation matrix dimensions
        oriDims = mxGetDimensions(mxOri);
        if (oriDims[0] != IM_NDIMS || oriDims[1] != IM_NDIMS)
                return SIFT3D_FAILURE;

        // Allocate space in the keypoint store
        if (resize_Keypoint_store(kp, (size_t) numKp))
                return SIFT3D_FAILURE;

        // Copy the data
        for (i = 0; i < (int) numKp; i++) {

                double *coordsData;

                Keypoint *const key = kp->buf + i;
                const mwIndex idx = (mwIndex) i;

                // Get the matrices
                if ((mxCoords = mxGetFieldByNumber(mx, idx, coordsNum)) == 
                                NULL ||
                        (mxScale = mxGetFieldByNumber(mx, idx, scaleNum)) == 
                                NULL ||
                        (mxOri = mxGetFieldByNumber(mx, idx, oriNum)) == 
                                NULL ||
                        (mxOctave = mxGetFieldByNumber(mx, idx, octaveNum)) == 
                                NULL ||
                        (mxLevel = mxGetFieldByNumber(mx, idx, levelNum)) == 
                                NULL)
                        return SIFT3D_FAILURE;


                // Copy the scalars
                key->sd = mxGetScalar(mxScale);
                key->o = (int) mxGetScalar(mxOctave);
                key->s = (int) mxGetScalar(mxLevel);

                // Get the coordinate data
                if ((coordsData = mxGetData(mxCoords)) == NULL)
                        return SIFT3D_FAILURE;

                // Copy the coordinate vector
                key->xd = coordsData[0];
                key->yd = coordsData[1];
                key->zd = coordsData[2];

		// Initialize the orientation matrix
		if (init_Mat_rm_p(&key->R, key->r_data, IM_NDIMS, IM_NDIMS, 
                        FLOAT, SIFT3D_FALSE))
			return SIFT3D_FAILURE;

                // Copy the orientation matrix
                if (mx2mat(mxOri, &key->R))
                        return SIFT3D_FAILURE;

                // Compute the remaining fields 
                key->xi = (int) key->xd;
                key->yi = (int) key->yd;
                key->zi = (int) key->zd;
                key->sd_rel = pow(2.0, -key->o) * key->sd;
        }

        return SIFT3D_SUCCESS;
}

/* Converts a SIFT3D_Descriptor_store to an mxArray.
 *
 * Output format:
 * [x y z el0 el1 ... el767] (see sift.c:SIFT3D_Descriptor_store_to_Mat_rm)
 *
 * Returns a pointer to the array, or NULL on failure. */
mxArray *desc2mx(const SIFT3D_Descriptor_store *const desc) {

        mxArray *mx;
        Mat_rm mat;

        // Initialize intermediates
        if (init_Mat_rm(&mat, 0, 0, FLOAT, SIFT3D_FALSE))
                return NULL;

        // Convert desc to a matrix
        if (SIFT3D_Descriptor_store_to_Mat_rm(desc, &mat))
                goto desc2mx_quit;

        // Convert the matrix to an mxArray
        if ((mx = mat2mx(&mat)) == NULL)
                goto desc2mx_quit;

        // Clean up
        cleanup_Mat_rm(&mat);

        return mx;

desc2mx_quit:
        cleanup_Mat_rm(&mat);
        return NULL;
} 

/* Wrapper for SIFT3D_detect_keypoints. */
int mex_SIFT3D_detect_keypoints(const Image *const im, 
        Keypoint_store *const kp) {
        return SIFT3D_detect_keypoints(&sift3d, im, kp);
}

/* Wrapper for SIFT3D_extract_descriptors. */
int mex_SIFT3D_extract_descriptors(const void *const im, 
        const Keypoint_store *const kp, SIFT3D_Descriptor_store *const desc, 
        const int useGpyr) {
        return SIFT3D_extract_descriptors(&sift3d, im, kp, desc, useGpyr);
}

/* Wrapper to get the Gaussian pyramid from the SIFT3D struct. */
Pyramid *mexGetGpyr(void) {
        return &sift3d.gpyr;
}
