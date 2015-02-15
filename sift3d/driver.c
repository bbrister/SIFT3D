/* driver.c
 * ----------------------------------------------------------------
 * Rice MRI Team
 * ----------------------------------------------------------------
 * This file contains the CLI to run SIFT on a volumetric image. It 
 * calls the routines in sift.c.
 *-----------------------------------------------------------------
 * Created: Blaine Rister 12/26/2013
 * Last updated: Blaine Rister 2/14/2013
 */

#include "stdio.h"
#include "types.h"
#include "sift.h"
#include "macros.h"
#include "imutil.h"

#define NN_THRESH 0.8

/* CLI for 3D SIFT. 
 * Format: sift3d [reference][source](arguments)
 * Arguments: 
 * --first_octave	  - the first octave 
 *						(Use -1000 for default)
 * --peak_thresh	  - threshold on DoG extrema magnitude
 *						(Use -1.0 for default)
 * --edge_thresh	  - threshold on edge score 
 *						(Use -1.0 for default)
 * --num_octaves	  - total number of octaves (default: process
 *						until one dimension is less than 8, use 
 *						-1 for default)
 * --num_kp_levels	  - number of levels per octave for keypoint
 *						candidates (Use -1 for default)
 * --sigma_n		  - base level of blurring assumed in data
 *						(default: Use -1.0 for default)
 */

int main(int argc, char *argv[]) {

	Image src, ref;
	SIFT3D_Detector detector;
	SIFT3D_Extractor extractor;
	Keypoint_store kp_src, kp_ref;
	SIFT3D_Descriptor_store desc_src, desc_ref;
	Mat_rm match_src, match_ref;
	int *matches;
	char *src_path, *ref_path;

	// FIXME: The CLI does not accept SIFT parameters
	if (argc < 3) {
		puts("Usage: sift3d [im1] [im2] \n");
		return 1;
	}

	src_path = argv[1];
	ref_path = argv[2];

	// Initialize data 
	init_Keypoint_store(&kp_src); 
	init_Keypoint_store(&kp_ref);
	init_SIFT3D_Descriptor_store(&desc_src); 
	init_SIFT3D_Descriptor_store(&desc_ref);
	init_im(&src);
	init_im(&ref);
	init_Mat_rm(&match_src, 0, 0, DOUBLE, FALSE);
	init_Mat_rm(&match_ref, 0, 0, DOUBLE, FALSE);
	matches = NULL;

	// Initialize the SIFT detector
	if (init_SIFT3D_Detector(&detector, 0))
		err_exit("init sift detector");

	// Initialize the SIFT descriptor extractor 
	if (init_SIFT3D_Extractor(&extractor))
		err_exit("init sift extractor");

	// Load the images
	if (read_nii(src_path, &src) ||
		read_nii(ref_path, &ref))
		err_exit("load images");

	// Extract features from reference image
	if (SIFT3D_detect_keypoints(&detector, &ref, &kp_ref))
		err_exit("detect reference keypoints");
	if (SIFT3D_extract_descriptors(&extractor, &detector.gpyr, &kp_ref,
		&desc_ref, TRUE))
		err_exit("extract reference descriptors");


	// Save intermediate data
#ifndef NDEBUG
	if (write_nii("../debug/src_im.nii", &src) ||
		write_nii("../debug/ref_im.nii", &ref))
		err_exit("write input images");
	if (write_Keypoint_store(KP_REF_PATH, &kp_ref))
		err_exit("write keypoints");
	printf("%u Reference keypoints written to %s\n",
		   (unsigned int) kp_ref.slab.num, KP_REF_PATH);
#endif
#ifdef VERBOSE
	if (write_pyramid(GPYR_REF_PATH, &detector.gpyr))
		fprintf(stderr, "Failed to write reference GSS pyramid to"
						" path %s \n", GPYR_REF_PATH);
		printf("Reference GSS pyramid written to %s \n", GPYR_REF_PATH);
	if (write_pyramid(DOG_REF_PATH, &detector.dog))
		fprintf(stderr, "Failed to write reference DoG pyramid to"
						" path %s \n", DOG_REF_PATH);
		printf("Reference DoG pyramid written to %s \n", DOG_REF_PATH);
#endif

	if (SIFT3D_detect_keypoints(&detector, &src, &kp_src))
		err_exit("detect source keypoints");
	if (SIFT3D_extract_descriptors(&extractor, &detector.gpyr, &kp_src, 
								   &desc_src, TRUE))
		err_exit("extract source descriptors");

#ifdef VERBOSE
	// More intermediate data
	if (write_pyramid(GPYR_SRC_PATH, &detector.gpyr))
		fprintf(stderr, "Failed to write source GSS pyramid to"
						" path %s \n", GPYR_SRC_PATH);
		printf("Source GSS pyramid written to %s \n", GPYR_SRC_PATH);
	if (write_pyramid(DOG_SRC_PATH, &detector.dog))
		fprintf(stderr, "Failed to write source DoG pyramid to"
						" path %s \n", DOG_SRC_PATH);
		printf("Source DoG pyramid written to %s \n", DOG_SRC_PATH);
#endif
#ifndef NDEBUG
	if (write_Keypoint_store(KP_SRC_PATH, &kp_src))
		err_exit("write keypoints");
	printf("%u Source keypoints written to %s\n",
		   (unsigned int) kp_src.slab.num, KP_SRC_PATH);
#endif

	// Match features and output results
	if (SIFT3D_nn_match_fb(&desc_src, &desc_ref, NN_THRESH, &matches)) 
		err_exit("match keypoints\n");

	if (SIFT3D_matches_to_Mat_rm(&desc_src, &desc_ref, matches,
				     &match_src, &match_ref))
		err_exit("extract coordinate matrices \n");

	if (write_Mat_rm(MATCH_SRC_PATH, &match_src) ||	
		write_Mat_rm(MATCH_REF_PATH, &match_ref))
		err_exit("write matches to file");
 	
	printf("%u matched features written to %s and %s \n", 
		   (unsigned int) match_src.num_rows, MATCH_SRC_PATH, 
		   MATCH_REF_PATH);

	return 0;
}
