/* sift.h
 * ----------------------------------------------------------------
 * Internal header for sift.c
 * ----------------------------------------------------------------
 */

#include "types.h"

#ifndef _SIFT_H
#define _SIFT_H

void init_Keypoint_store(Keypoint_store *kp);

void init_SIFT3D_Descriptor_store(SIFT3D_Descriptor_store *desc);

int set_first_octave_SIFT3D(SIFT3D *const sift3d, 
                                const int first_octave);

int set_peak_thresh_SIFT3D(SIFT3D *const sift3d,
                                const double peak_thresh);

int set_corner_thresh_SIFT3D(SIFT3D *const sift3d,
                                const double corner_thresh);

int set_num_octaves_SIFT3D(SIFT3D *const sift3d,
                                const unsigned int num_octaves);

int set_num_kp_levels_SIFT3D(SIFT3D *const sift3d,
                                const unsigned int num_kp_levels);

int set_sigma_n_SIFT3D(SIFT3D *const sift3d,
                                const double sigma_n);

int set_sigma0_SIFT3D(SIFT3D *const sift3d,
                                const double sigma_n);

int init_SIFT3D(SIFT3D *sift3d);

void print_opts_SIFT3D(void);

int parse_args_SIFT3D(SIFT3D *const sift3d,
        const int argc, char *const *argv, int *optind_ret, 
        const int check_err);

int SIFT3D_detect_keypoints(SIFT3D *const sift3d, Image *const im, 
			    Keypoint_store *const kp);

int SIFT3D_extract_descriptors(SIFT3D *sift3d, void *im,
 	Keypoint_store *kp, SIFT3D_Descriptor_store *desc, int use_gpyr);

int SIFT3D_extract_dense_descriptors(SIFT3D *const sift3d, 
        const Image *const in, Image *const desc);

int SIFT3D_nn_match(const SIFT3D_Descriptor_store *const d1,
		    const SIFT3D_Descriptor_store *const d2,
		    const float nn_thresh, int **const matches);

int SIFT3D_nn_match_fb(const SIFT3D_Descriptor_store *const d1,
		       const SIFT3D_Descriptor_store *const d2,
		       const float nn_thresh, int **const matches);

int Keypoint_store_to_Mat_rm(const Keypoint_store *const kp, Mat_rm *const mat);

int SIFT3D_Descriptor_store_to_Mat_rm(const SIFT3D_Descriptor_store *const store, 
				      Mat_rm *const mat);

int SIFT3D_matches_to_Mat_rm(SIFT3D_Descriptor_store *d1,
			     SIFT3D_Descriptor_store *d2,
			     const int *const matches,
			     Mat_rm *const match1, 
			     Mat_rm *const match2);

int draw_matches(const Image *const src, const Image *const ref,
		 const Mat_rm *const match_src, const Mat_rm *const match_ref,
		 Image *const background, Image *const overlay);

#endif
