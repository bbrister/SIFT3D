/* -----------------------------------------------------------------------------
 * mexutil.h
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Public header for SIFT3D mex utility library.
 * -----------------------------------------------------------------------------
 */

#include <uchar.h>
#include "mex.h"
#include "imtypes.h"

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

int mx2desc(const mxArray *const mx, SIFT3D_Descriptor_store *const desc);

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

int mex_set_opts_SIFT3D(const mxArray *const mx);

int mex_set_opts_Reg_SIFT3D(const mxArray *const mx);

int mex_set_nn_thresh_Reg_SIFT3D(const double nn_thresh);

double mex_get_nn_thresh_Reg_SIFT3D(void);

int mex_register_SIFT3D_resample(const Image *const src, 
        const Image *const ref, const interp_type interp, void *const tform);

int mex_set_src_Reg_SIFT3D(const Image *const src);

int mex_set_ref_Reg_SIFT3D(const Image *const ref);

int mex_register_SIFT3D(void *const tform);

int mex_get_matches_Reg_SIFT3D(Mat_rm *const match_src, 
        Mat_rm *const match_ref);

#ifdef __cplusplus
}
#endif

#endif
