/* reg.c
 * ----------------------------------------------------------------
 * This files contains routines for image registraion. 
 * ----------------------------------------------------------------
 */

#include "types.h"

#ifndef _REG_H
#define _REG_H

/* Internal data for the SIFT3D + RANSAC registration process */
typedef struct _Reg_SIFT3D {

        SIFT3D sift3d;
        Ransac ran;
        Image src, ref;
        Keypoint_store kp_src, kp_ref;
        SIFT3D_Descriptor_store desc_src, desc_ref;
        Mat_rm match_src, match_ref;
        int *matches;
        double nn_thresh;
        int verbose;

} Reg_SIFT3D;

int init_Reg_SIFT3D(Reg_SIFT3D *const reg);

int register_SIFT3D(Reg_SIFT3D *const reg, const tform_type mode, 
        void *const tform);

int set_src_Reg_SIFT3D(Reg_SIFT3D *const reg, Image *const src);

int set_ref_Reg_SIFT3D(Reg_SIFT3D *const reg, Image *const ref);

#endif
