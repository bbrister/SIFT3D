/* -----------------------------------------------------------------------------
 * mexutil.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Mex utility library for SIFT3D.
 * -----------------------------------------------------------------------------
 */

/* Standard headers */
#include <stddef.h>
#include <math.h>

/* SIFT3D headers */
#include "imtypes.h"
#include "immacros.h"
#include "imutil.h"
#include "sift.h"
#include "reg.h"
#include "mexutil.h"

/* Matlab headers */
#include <uchar.h>
#include "mex.h"
#include "matrix.h"

/* The number of dimensions in mxArrays representing images */
#define MX_IM_NDIMS (IM_NDIMS + 1) 

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

/* Global state */
Reg_SIFT3D reg;

/* Error message tag */
const char *tag = "sift3D";

/* Static helper functions */
static void init(void) __attribute__((constructor));
static void fini(void) __attribute__((destructor));

/* Library initialization */
static void init(void) {
        if (init_Reg_SIFT3D(&reg))
                err_msgu("main:initSift", "Failed to initialize SIFT3D");
}

/* Library cleanup */
static void fini(void) {
        cleanup_Reg_SIFT3D(&reg);
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

/* Returns the index in an mxArray of the voxel at the coordinates (x, y, z)
 * and channel c */
mwIndex mxImGetIdx(const mxArray *const mx, const int x, const int y, 
        const int z, const int c) {

        const mwIndex subs[] = {x, y, z, c};
        const mwSize nSubs = sizeof(subs) / sizeof(mwIndex);
        assert(nSubs == MX_IM_NDIMS);

        return mxCalcSingleSubscript(mx, nSubs, subs);
}

/* Convert an Image to an mxArray. The output will have IM_NDIMS + 1 dimensions
 * and type double. The final dimension denotes the image channels. 
 *
 * Returns a pointer to the array, or NULL if an error has occurred. */
mxArray *im2mx(const Image *const im) {
    
        mwSize dims[MX_IM_NDIMS]; 
        mxArray *mx; 
        double *mxData; 
        int i, x, y, z, c;

        // Initialize the dimensions
        for (i = 0; i < IM_NDIMS; i++) {
                dims[i] = SIFT3D_IM_GET_DIMS(im)[i];
        }
        dims[IM_NDIMS] = im->nc;

        // Create an array
        if ((mx = mxCreateNumericArray(MX_IM_NDIMS, dims, mxDOUBLE_CLASS, 
                mxREAL)) == NULL)
                return NULL;

        // Get the data
        if ((mxData = mxGetData(mx)) == NULL)
                return NULL;

        // Transpose and copy the data
        SIFT3D_IM_LOOP_START_C(im, x, y, z, c)

                const mwIndex idx = mxImGetIdx(mx, x, y, z, c); 
                mxData[idx] = (double) SIFT3D_IM_GET_VOX(im, x, y, z, c);

        SIFT3D_IM_LOOP_END_C

        return mx;
}

/* Convert an mxArray to an Image. mx must have at most IM_NDIMS + 1 
 * dimensions and be of type single (float). */
int mx2im(const mxArray *const mx, Image *const im) {

        const mwSize *mxDims;
        float *mxData;
        mwIndex mxNDims;
        int i, x, y, z, c, mxNSpaceDims;

        // Verify inputs
        mxNDims = (int) mxGetNumberOfDimensions(mx);
	if (mxNDims > MX_IM_NDIMS || !mxIsSingle(mx) || mxIsComplex(mx))
                return SIFT3D_FAILURE;

        // Copy the spatial dimensions
        mxNSpaceDims = SIFT3D_MIN((int) mxNDims, MX_IM_NDIMS - 1);
        mxDims = mxGetDimensions(mx);
        for (i = 0; i < mxNSpaceDims; i++) {
                SIFT3D_IM_GET_DIMS(im)[i] = (int) mxDims[i];
        }

        // Pad the unfilled dimensions with 1
        for (i = mxNSpaceDims; i < MX_IM_NDIMS - 1; i++) {
                SIFT3D_IM_GET_DIMS(im)[i] = 1; 
        } 

        // Copy the number of channels, defaulting to 1
        im->nc = mxNDims == MX_IM_NDIMS ? (int) mxDims[MX_IM_NDIMS - 1] : 1;

        // Resize the output                
        im_default_stride(im);
        if (im_resize(im))
                return SIFT3D_FAILURE;

        // Get the data
        if ((mxData = mxGetData(mx)) == NULL)
                return SIFT3D_FAILURE;

        // Transpose and copy the data
        SIFT3D_IM_LOOP_START_C(im, x, y, z, c)

                mwIndex idx;

                idx = mxImGetIdx(mx, x, y, z, c);

                SIFT3D_IM_GET_VOX(im, x, y, z, c) = mxData[idx];

        SIFT3D_IM_LOOP_END_C

        return SIFT3D_SUCCESS;
}

/* Returns an mxArray representing the units of image im.
 *
 * Parameters:
 * Return: the array, or NULL on failure. */
mxArray *units2mx(const Image *const im) {

        mxArray *mx;
        double *mxData;
        int i;

        // Create an array
        if ((mx = mxCreateDoubleMatrix(IM_NDIMS, 1, mxREAL)) == NULL)
                return NULL;

        // Get the data
        if ((mxData = mxGetData(mx)) == NULL)
                goto units2mx_quit;

        // Copy the data from im
        for (i = 0; i < IM_NDIMS; i++) {

                mwIndex idx;

                const mwIndex subs[] = {i, 0};
                const mwIndex nSubs = sizeof(subs) / sizeof(mwIndex);

                idx = mxCalcSingleSubscript(mx, nSubs, subs);
                mxData[idx] = SIFT3D_IM_GET_UNITS(im)[i];
        }

        return mx;

units2mx_quit:
        mxDestroyArray(mx);
        return NULL;
}

/* Set the units of an image to the data stored in an mxArray,
 *
 * Parameters:
 *  -mx: An array of IM_NDIMS dimensions, type double.
 *  -im: The Image struct to which the units are written.
 *
 * Return: SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int mx2units(const mxArray *const mx, Image *const im) {

        const mwSize *mxDims;
        double *mxData;
        mwIndex mxNDims;
        int i;

        // Verify inputs
        mxNDims = (int) mxGetNumberOfDimensions(mx);
	if (mxNDims != 2 || !isDouble(mx))
                return SIFT3D_FAILURE;

        // Verify the dimensions
        mxDims = mxGetDimensions(mx);
        if (mxDims[0] != IM_NDIMS || mxDims[1] != 1)
                return SIFT3D_FAILURE; 

        // Get the data
        if ((mxData = mxGetData(mx)) == NULL)
                return SIFT3D_FAILURE;

        // Copy the data to im
        for (i = 0; i < IM_NDIMS; i++) {

                mwIndex idx;

                const mwIndex subs[] = {i, 0};
                const mwIndex nSubs = sizeof(subs) / sizeof(mwIndex);

                idx = mxCalcSingleSubscript(mx, nSubs, subs);
                SIFT3D_IM_GET_UNITS(im)[i] = mxData[idx];
        }

        return SIFT3D_SUCCESS;
}

/* Wrapper around mx2im and mx2units. */
int mx2imWithUnits(const mxArray *const data, const mxArray *const units,
        Image *const im) {
        return mx2im(data, im) || mx2units(units, im);
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
        SIFT3D_MAT_RM_TYPE_MACRO(mat, mat2mx_quit, TRANSPOSE_AND_COPY);

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
        if (!isDouble(mx)) {
                SIFT3D_ERR("mx2mat: mx must have type double");
                return SIFT3D_FAILURE;
        }

        if (mxGetNumberOfDimensions(mx) != 2) {
                SIFT3D_ERR("mx2mat: mx must be a matrix");
                return SIFT3D_FAILURE;
        }

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
        SIFT3D_MAT_RM_TYPE_MACRO(mat, mx2mat_quit, COPY_DATA);

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
                        SIFT3D_FLOAT, SIFT3D_FALSE))
			return SIFT3D_FAILURE;

                // Copy the orientation matrix
                if (mx2mat(mxOri, &key->R))
                        return SIFT3D_FAILURE;
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
        if (init_Mat_rm(&mat, 0, 0, SIFT3D_FLOAT, SIFT3D_FALSE))
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

/* Converts a pair of mxArrays to a SIFT3D_Descriptor_store. */
int mx2desc(const mxArray *const mx, SIFT3D_Descriptor_store *const desc) {

        Mat_rm mat;

        // Initialize intermediates
        if (init_Mat_rm(&mat, 0, 0, SIFT3D_FLOAT, SIFT3D_FALSE))
                return SIFT3D_FAILURE;

        // Convert the mxArrays to a matrix
        if (mx2mat(mx, &mat)) {
                SIFT3D_ERR("mx2desc: Failed to convert mx to mat");
                goto mx2desc_quit;
        }

        // Convert the matrix to a SIFT3D_Descriptor_store
        if (Mat_rm_to_SIFT3D_Descriptor_store(&mat, desc)) {
                SIFT3D_ERR("mx2desc: Failed to convert mat to desc");
                goto mx2desc_quit;
        }

        // Clean up
        cleanup_Mat_rm(&mat);

        return SIFT3D_SUCCESS;

mx2desc_quit:
        cleanup_Mat_rm(&mat);
        return SIFT3D_FAILURE;
}

/* Convert an array of doubles to an mxArray. 
 *
 * Parameters:
 *   array: The array.
 *   num: The number of elements.
 *
 * Returns a pointer to the mxArray, or NULL on failure. */
mxArray *array2mx(const double *const array, const size_t len) {

        mxArray *mx;
        double *mxData; 
        size_t i;

        const mwSize dims[] = {len}; 
        const int ndim = 1; 

        // Create the mxArray
        if ((mx = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, 
                mxREAL)) == NULL)
                return NULL;

        // Get a pointer to the array's data
        if ((mxData = mxGetData(mx)) == NULL)
                return NULL;

        // Copy the data
        for (i = 0; i < len; i++) {
                mxData[i] = array[i];
        }

        return mx;
}

/* Wrapper for SIFT3D_detect_keypoints. */
int mex_SIFT3D_detect_keypoints(const Image *const im, 
        Keypoint_store *const kp) {
        return SIFT3D_detect_keypoints(&reg.sift3d, im, kp);
}

/* Wrapper for SIFT3D_assign_orientations. */
int mex_SIFT3D_assign_orientations(const Image *const im, 
        Keypoint_store *const kp, double **const conf) {
        return SIFT3D_assign_orientations(&reg.sift3d, im, kp, conf);
}

/* Wrapper for SIFT3D_extract_descriptors. */
int mex_SIFT3D_extract_descriptors(const Keypoint_store *const kp, 
        SIFT3D_Descriptor_store *const desc) {
        return SIFT3D_extract_descriptors(&reg.sift3d, kp, desc);
}

/* Wrapper for SIFT3D_extract_raw_descriptors. */
int mex_SIFT3D_extract_raw_descriptors(const Image *const im, 
        const Keypoint_store *const kp, SIFT3D_Descriptor_store *const desc) {
        return SIFT3D_extract_raw_descriptors(&reg.sift3d, im, kp, desc);
}

/* Wrapper for SIFT3D_have_gpyr. */
int mexHaveGpyr(void) {
        return SIFT3D_have_gpyr(&reg.sift3d);
}

/* Wrapper to set the options for the SIFT3D struct. The argument mx shall be
 * a struct with the following fields:
 *      -peakThresh
 *      -cornerThresh
 *      -numKpLevels
 *      -sigmaN
 *      -sigma0
 *
 * Any other fields are ignored. Empty fields are replaced with the defaults. 
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int mex_set_opts_SIFT3D(const mxArray *const mx) {

        SIFT3D sift3d;
        const mxArray *mxFirstOctave, *mxPeakThresh, *mxCornerThresh, 
                *mxNumKpLevels, *mxSigmaN, *mxSigma0;

        // Verify inputs
        if (mxIsEmpty(mx) || !mxIsStruct(mx))
                return SIFT3D_FAILURE;

        // Get the option arrays
        if ((mxPeakThresh = mxGetField(mx, 0, "peakThresh")) == NULL ||
                (mxCornerThresh = mxGetField(mx, 0, "cornerThresh")) == NULL ||
                (mxNumKpLevels = mxGetField(mx, 0, "numKpLevels")) == NULL ||
                (mxSigmaN = mxGetField(mx, 0, "sigmaN")) == NULL ||
                (mxSigma0 = mxGetField(mx, 0, "sigma0")) == NULL)
                return SIFT3D_FAILURE;

        // Initialize intermediates
        if (init_SIFT3D(&sift3d))
                return SIFT3D_FAILURE;

        // Set the non-empty options in our new SIFT3D struct
        if (!mxIsEmpty(mxPeakThresh) && 
                set_peak_thresh_SIFT3D(&sift3d, mxGetScalar(mxPeakThresh)))
                goto mex_set_opts_SIFT3D_quit;
        if (!mxIsEmpty(mxCornerThresh) &&
                set_corner_thresh_SIFT3D(&sift3d, mxGetScalar(mxCornerThresh)))
                goto mex_set_opts_SIFT3D_quit;
        if (!mxIsEmpty(mxNumKpLevels) &&
                set_num_kp_levels_SIFT3D(&sift3d, 
                (int) mxGetScalar(mxNumKpLevels)))
                goto mex_set_opts_SIFT3D_quit;
        if (!mxIsEmpty(mxSigmaN) &&
                set_sigma_n_SIFT3D(&sift3d, mxGetScalar(mxSigmaN)))
                goto mex_set_opts_SIFT3D_quit;
        if (!mxIsEmpty(mxSigma0) &&
                set_sigma0_SIFT3D(&sift3d, mxGetScalar(mxSigma0)))
                goto mex_set_opts_SIFT3D_quit;

        // Set the options in the Reg_SIFT3D struct
        if (set_SIFT3D_Reg_SIFT3D(&reg, &sift3d))
                goto mex_set_opts_SIFT3D_quit;

        // Clean up
        cleanup_SIFT3D(&sift3d);

        return SIFT3D_SUCCESS;

mex_set_opts_SIFT3D_quit:
        cleanup_SIFT3D(&sift3d);
        return SIFT3D_FAILURE;
}

/* Wrapper to set the options for the Reg_SIFT3D struct. The argument mx shall
 * be a struct with the following fields:
 *   -numIter
 *   -errThresh
 *   -nnThresh
 *
 * Any other fields are ignored. Empty fields are replaced with the defaults.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int mex_set_opts_Reg_SIFT3D(const mxArray *const mx) {

        Ransac ran;
        const mxArray *mxNumIter, *mxErrThresh, *mxNnThresh;
        double nn_thresh;

        // Verify inputs
        if (mxIsEmpty(mx) || !mxIsStruct(mx))
                return SIFT3D_FAILURE;
        
        // Initialize options to the defaults
        nn_thresh = SIFT3D_nn_thresh_default;
        init_Ransac(&ran);

        // Get the option arrays
        if ((mxNumIter = mxGetField(mx, 0, "numIter")) == NULL ||
                (mxErrThresh = mxGetField(mx, 0, "errThresh")) == NULL ||
                (mxNnThresh = mxGetField(mx, 0, "nnThresh")) == NULL)
                return SIFT3D_FAILURE;

        // Get the non-empty options 
        if (!mxIsEmpty(mxNumIter) &&
                set_num_iter_Ransac(&ran, (int) mxGetScalar(mxNumIter))) {
                return SIFT3D_FAILURE;
        } 
        if (!mxIsEmpty(mxErrThresh) && 
                set_err_thresh_Ransac(&ran, mxGetScalar(mxErrThresh))) {
                return SIFT3D_FAILURE;
        } 
        if (!mxIsEmpty(mxNnThresh)) {
                nn_thresh = mxGetScalar(mxNnThresh);
        }

        // Set the options
        return set_Ransac_Reg_SIFT3D(&reg, &ran) ||
            set_nn_thresh_Reg_SIFT3D(&reg, nn_thresh);
}

/* Get the current value of nn_thresh from the Reg_SIFT3D struct. */
double mex_get_nn_thresh_Reg_SIFT3D(void) {
        return reg.nn_thresh;
}

/* Wrapper for register_SIFT3D_resample */
int mex_register_SIFT3D_resample(const Image *const src, 
        const Image *const ref, const interp_type interp, void *const tform) {
        return register_SIFT3D_resample(&reg, src, ref, interp, tform); 
}

/* Wrapper for set_src_Reg_SIFT3D */
int mex_set_src_Reg_SIFT3D(const Image *const src) {
        return set_src_Reg_SIFT3D(&reg, src); 
}

/* Wrapper for set_ref_Reg_SIFT3D */
int mex_set_ref_Reg_SIFT3D(const Image *const ref) {
        return set_ref_Reg_SIFT3D(&reg, ref); 
}

/* Wrapper for register_SIFT3D */
int mex_register_SIFT3D(void *const tform) {
        return register_SIFT3D(&reg, tform);
}

/* Wrapper for get_matches_Reg_SIFT3D */
int mex_get_matches_Reg_SIFT3D(Mat_rm *const match_src, 
        Mat_rm *const match_ref) {
        return get_matches_Reg_SIFT3D(&reg, match_src, match_ref);
}
