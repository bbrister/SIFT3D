/* -----------------------------------------------------------------------------
 * reg.h
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Public header file for reg.c.
 * -----------------------------------------------------------------------------
 */

#include "imtypes.h"

#ifndef _REG_H
#define _REG_H

#ifdef __cplusplus
extern "C" {
#endif

/* Parameters */
const extern double SIFT3D_nn_thresh_default; // Default matching threshold

/* Internal data for the SIFT3D + RANSAC registration process */
typedef struct _Reg_SIFT3D {

        double src_units[IM_NDIMS], ref_units[IM_NDIMS];
        SIFT3D sift3d;
        Ransac ran;
        SIFT3D_Descriptor_store desc_src, desc_ref;
        Mat_rm match_src, match_ref;
        double nn_thresh;
        int verbose;

} Reg_SIFT3D;

int init_Reg_SIFT3D(Reg_SIFT3D *const reg);

void cleanup_Reg_SIFT3D(Reg_SIFT3D *const reg);

int register_SIFT3D(Reg_SIFT3D *const reg, void *const tform);

int register_SIFT3D_resample(Reg_SIFT3D *const reg, const Image *const src,
	const Image *const ref, const interp_type interp, void *const tform);

int set_nn_thresh_Reg_SIFT3D(Reg_SIFT3D *const reg, const double nn_thresh);

int set_Ransac_Reg_SIFT3D(Reg_SIFT3D *const reg, const Ransac *const ran);

int set_SIFT3D_Reg_SIFT3D(Reg_SIFT3D *const reg, const SIFT3D *const sift3d);

int set_src_Reg_SIFT3D(Reg_SIFT3D *const reg, const Image *const src);

int set_ref_Reg_SIFT3D(Reg_SIFT3D *const reg, const Image *const ref);

int get_matches_Reg_SIFT3D(const Reg_SIFT3D *const reg, Mat_rm *const match_src,
        Mat_rm *const match_ref);

#ifdef __cplusplus
}
#endif

#endif
