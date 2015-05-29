/* sift.h
 * ----------------------------------------------------------------
 * Rice MRI Team
 * ----------------------------------------------------------------
 * Internal header for sift.c
 *-----------------------------------------------------------------
 * Created: Blaine Rister 12/26/2013
 * Last updated: Blaine Rister 2/19/2014
 */

#include "types.h"

#ifndef _SIFT_H
#define _SIFT_H

/* Default SIFT detector parameters. These can be overriden by 
 * the appropriate initializer arguments. */
#define DEFAULT_FIRST_OCTAVE 0
#define DEFAULT_PEAK_THRESH 0.03
#define DEFAULT_CORNER_THRESH 0.1
#define DEFAULT_NUM_KP_LEVELS 3
#define DEFAULT_SIGMA_N 1.15 
#define DEFAULT_SIGMA0 1.6 

/* Internal parameters */
//#define USE_OPENCL				// Use GPU acceleration
#define PARABOLOID_INTERP			// Interpolate sub-bin orientations
//#define ORI_SOLID_ANGLE_WEIGHT	// Weight bins by solid angle
//#define MATCH_MAX_DIST 0.3 // Maximum distance between matching features 
						   // (ratio with volume diagonal)	
#define EIG_ORI 			// Use eigenvalue orientation algorithm

/* Safety checks */
#if defined ICOS_HIST && !defined EIG_ORI
#pragma error("Cannot use ICOS_HIST without EIG_ORI")
#endif

/* Debugging information */
#define DEBUG_ROOT "../debug/"
#define GPYR_SRC_PATH DEBUG_ROOT "gpyr_src/gpyr"
#define GPYR_REF_PATH DEBUG_ROOT "gpyr_ref/gpyr"
#define DOG_SRC_PATH DEBUG_ROOT "dog_src/dog"
#define DOG_REF_PATH DEBUG_ROOT "dog_ref/dog"
#define KP_SRC_PATH DEBUG_ROOT "kp_src.txt"
#define KP_REF_PATH DEBUG_ROOT "kp_ref.txt"
#define MATCH_SRC_PATH DEBUG_ROOT "match_src.txt"
#define MATCH_REF_PATH DEBUG_ROOT "match_ref.txt"

void init_Keypoint_store(Keypoint_store *kp);

void init_SIFT3D_Descriptor_store(SIFT3D_Descriptor_store *desc);

void set_first_octave_SIFT3D(SIFT3D *const sift3d, 
                                const int first_octave);

void set_peak_thresh_SIFT3D(SIFT3D *const sift3d,
                                const double peak_thresh);

void set_corner_thresh_SIFT3D(SIFT3D *const sift3d,
                                const double corner_thresh);

int set_num_octaves_SIFT3D(SIFT3D *const sift3d,
                                const unsigned int num_octaves);

int set_num_kp_levels_SIFT3D(SIFT3D *const sift3d,
                                const unsigned int num_kp_levels);

void set_sigma_n_SIFT3D(SIFT3D *const sift3d,
                                const double sigma_n);

void set_sigma0_SIFT3D(SIFT3D *const sift3d,
                                const double sigma_n);

int init_SIFT3D(SIFT3D *sift3d);

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
