/* -----------------------------------------------------------------------------
 * mexutil.h
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Public header for SIFT3D mex utility library.
 * -----------------------------------------------------------------------------
 */

#include "mex.h"
#include "types.h"

#ifndef _MEXUTIL_H
#define _MEXUTIL_H

#ifdef __cplusplus
extern "C" {
#endif

void err_msg(const char *name, const char *msg);

void err_msgu(const char *name, const char *msg);

int isDouble(const mxArray *const mx);

mwIndex mxImGetIdx(const mxArray *const mx, const int x, const int y, 
        const int z, const int c);

mxArray *im2mx(const Image *const im);

int mx2im(const mxArray *const mx, Image *const im);

mxArray *units2mx(const Image *const im);

int mx2units(const mxArray *const mx, Image *const im);

int mx2imWithUnits(const mxArray *const data, const mxArray *const units,
        Image *const im);

mxArray *mat2mx(const Mat_rm *const mat);

int mx2mat(const mxArray *const mx, Mat_rm *const mat);

mxArray *kp2mx(const Keypoint_store *const);

int mx2kp(const mxArray *const mx, Keypoint_store *const);

mxArray *desc2mx(const SIFT3D_Descriptor_store *const desc);

mxArray *array2mx(const double *const array, const size_t len);

int mex_SIFT3D_detect_keypoints(const Image *const im, 
        Keypoint_store *const kp);

int mex_SIFT3D_assign_orientations(const Image *const im, 
        Keypoint_store *const kp, double **const conf);

int mex_SIFT3D_extract_descriptors(const Keypoint_store *const kp, 
        SIFT3D_Descriptor_store *const desc);

int mex_SIFT3D_extract_raw_descriptors(const Image *const im, 
        const Keypoint_store *const kp, SIFT3D_Descriptor_store *const desc);

int mexHaveGpyr(void);

#ifdef __cplusplus
}
#endif

#endif
