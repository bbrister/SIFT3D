/* -----------------------------------------------------------------------------
 * sift.h
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Public header for sift.c
 * -----------------------------------------------------------------------------
 */

#include "imtypes.h"

#ifndef _SIFT_H
#define _SIFT_H

#ifdef __cplusplus
extern "C" {
#endif

void init_Keypoint_store(Keypoint_store *const kp);

int init_Keypoint(Keypoint *const key);

int resize_Keypoint_store(Keypoint_store *const kp, const size_t num);

int copy_Keypoint(const Keypoint *const src, Keypoint *const dst);

void cleanup_Keypoint_store(Keypoint_store *const kp);

void init_SIFT3D_Descriptor_store(SIFT3D_Descriptor_store *const desc);

void cleanup_SIFT3D_Descriptor_store(SIFT3D_Descriptor_store *const desc);

int set_peak_thresh_SIFT3D(SIFT3D *const sift3d,
                                const double peak_thresh);

int set_corner_thresh_SIFT3D(SIFT3D *const sift3d,
                                const double corner_thresh);

int set_num_kp_levels_SIFT3D(SIFT3D *const sift3d,
                                const unsigned int num_kp_levels);

int set_sigma_n_SIFT3D(SIFT3D *const sift3d,
                                const double sigma_n);

int set_sigma0_SIFT3D(SIFT3D *const sift3d,
                                const double sigma_n);

int init_SIFT3D(SIFT3D *sift3d);

int copy_SIFT3D(const SIFT3D *const src, SIFT3D *const dst);

void cleanup_SIFT3D(SIFT3D *const sift3d);

void print_opts_SIFT3D(void);

int parse_args_SIFT3D(SIFT3D *const sift3d,
        const int argc, char **argv, const int check_err);

int SIFT3D_assign_orientations(const SIFT3D *const sift3d, 
        const Image *const im, Keypoint_store *const kp, double **const conf);

int SIFT3D_detect_keypoints(SIFT3D *const sift3d, const Image *const im,
			    Keypoint_store *const kp);

int SIFT3D_have_gpyr(const SIFT3D *const sift3d);

int SIFT3D_extract_descriptors(SIFT3D *const sift3d, 
        const Keypoint_store *const kp, 
        SIFT3D_Descriptor_store *const desc);

int SIFT3D_extract_raw_descriptors(SIFT3D *const sift3d, 
        const Image *const im, const Keypoint_store *const kp, 
        SIFT3D_Descriptor_store *const desc);

int SIFT3D_extract_dense_descriptors(SIFT3D *const sift3d, 
        const Image *const in, Image *const desc);

int SIFT3D_nn_match(const SIFT3D_Descriptor_store *const d1,
		    const SIFT3D_Descriptor_store *const d2,
		    const float nn_thresh, int **const matches);

int Keypoint_store_to_Mat_rm(const Keypoint_store *const kp, Mat_rm *const mat);

int SIFT3D_Descriptor_coords_to_Mat_rm(
        const SIFT3D_Descriptor_store *const store, 
        Mat_rm *const mat);

int SIFT3D_Descriptor_store_to_Mat_rm(const SIFT3D_Descriptor_store *const store, 
				      Mat_rm *const mat);

int Mat_rm_to_SIFT3D_Descriptor_store(const Mat_rm *const mat, 
				      SIFT3D_Descriptor_store *const store);

int SIFT3D_matches_to_Mat_rm(SIFT3D_Descriptor_store *d1,
			     SIFT3D_Descriptor_store *d2,
			     const int *const matches,
			     Mat_rm *const match1, 
			     Mat_rm *const match2);

int draw_matches(const Image *const left, const Image *const right,
                 const Mat_rm *const keys_left, const Mat_rm *const keys_right,
		 const Mat_rm *const match_left, const Mat_rm *const match_right,
		 Image *const concat, Image *const keys, Image *const lines);

int write_Keypoint_store(const char *path, const Keypoint_store *const kp);

int write_SIFT3D_Descriptor_store(const char *path, 
        const SIFT3D_Descriptor_store *const desc);

#ifdef __cplusplus
}
#endif

#endif
