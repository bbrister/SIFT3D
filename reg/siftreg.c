/* siftreg.c
 * ----------------------------------------------------------------
 * Rice MRI Team
 * ----------------------------------------------------------------
 * Complete registration pipeline with synthetically generated 
 * source images.
 * ----------------------------------------------------------------
 * Created: Blaine Rister 2/17/2014
 * Last updated: Blaine Rister 11/18/2014
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <unistd.h>
#include "types.h"
#include "sift.h"
#include "macros.h"
#include "imutil.h"
#include "math.h"

/* Parameters */
#define INTERP LINEAR // Interpolation type for synthetic transform
#define NN_THRESH_DEFAULT 0.8
#define MIN_INLIERS 0.001
#define ERR_THRESH 2.0
#define NUM_ITER 500

/* Debugging file paths */
#define TFORM_IN_PATH DEBUG_ROOT "affine_input.m"
#define TFORM_REG_PATH DEBUG_ROOT "affine_inferred.m"
#define IM_OUT_PATH DEBUG_ROOT "registered.nii"
#define GRID_PATH DEBUG_ROOT "grid.nii"
#define KP_SRC_IM_PATH DEBUG_ROOT "kp_src.nii"
#define KP_REF_IM_PATH DEBUG_ROOT "kp_ref.nii"
#define BG_PATH DEBUG_ROOT "background.nii"
#define OVERLAY_PATH DEBUG_ROOT "overlay.nii"

int main(int argc, char *argv[]) {

	Image src, ref, srcp, refp, srcp_reg;
	SIFT3D_Detector detector;
	SIFT3D_Extractor extractor;
	Keypoint_store kp_src, kp_ref;
	SIFT3D_Descriptor_store desc_src, desc_ref;
	Mat_rm match_src, match_ref, A;
	Ransac ran;
	Affine aff_syn, aff_reg;
	int *matches;
	char *ref_path, *src_path;
	double nn_thresh;
	struct timespec reg_start, reg_end, diff;
	long elapsed_ms;
	int i, c, syn_mode, nx, ny, nz;

	// Default parameters
	nn_thresh = NN_THRESH_DEFAULT;
	syn_mode = 0;

	// Initialize the SIFT detector
	if (init_SIFT3D_Detector(&detector))
		err_exit("init sift detector\n");

	// Initialize the SIFT descriptor extractor 
	if (init_SIFT3D_Extractor(&extractor))
		err_exit("init sift extractor\n");

        // Parse the SIFT3D options and increment the argument list
        parse_args_SIFT3D_Detector(&detector, argc, argv, &optind, 0);

	// Parse the options	
	while ((c = getopt(argc, argv, "+sn:")) != -1)
		switch (c) {
		case 'n':
			nn_thresh = atof(optarg);
			break;
		case 's':
			syn_mode = 1;
			break;
	        case '?':
			if (optopt == 'n')
			  fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else if (isprint(optopt))
			  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
			else
			  fprintf (stderr,
				   "Unknown option character `\\x%x'.\n",
				   optopt);
			return 1;
                        break;
	      	default:
			return 1;
	}

	// Parse the required arguments
	switch (argc - optind)  {
	case 1:
		puts("Registering synthetic data \n"); 
		ref_path = argv[optind];
		src_path = NULL;
		syn_mode = 1;
		break;
	case 2:
		puts("Registering real data \n"); 
		ref_path = argv[optind];
		src_path = argv[optind + 1];
		break;
      	default:
		puts("Usage: siftreg [reference image] [source image] \n" 
			 "Only use im1 for synthetic \n");
		return 1;
	}
	printf("nn_thresh: %f \n syn_mode: %d \n", nn_thresh, syn_mode);

	// Initialize data
	matches = NULL;
	init_Keypoint_store(&kp_src); 
	init_Keypoint_store(&kp_ref);
	init_SIFT3D_Descriptor_store(&desc_src); 
	init_SIFT3D_Descriptor_store(&desc_ref);
	init_im(&src);
	init_im(&srcp);
	init_im(&ref);
	init_im(&refp);
	init_im(&srcp_reg);
	init_Affine(&aff_syn, 3);
	init_Affine(&aff_reg, 3);
	init_Ransac(&ran);
	if (init_Mat_rm(&match_src, 0, 0, DOUBLE, FALSE) ||
		init_Mat_rm(&match_ref, 0, 0, DOUBLE, FALSE) ||
		init_Mat_rm(&A, 3, 4, DOUBLE, TRUE))
		err_exit("init data");

        // Set the RANSAC parameters
        set_min_inliers_Ransac(&ran, MIN_INLIERS);
        set_err_thresh_Ransac(&ran, ERR_THRESH);
        set_num_iter_Ransac(&ran, NUM_ITER);

	// Load the reference image
	if (read_nii(ref_path, &ref))
		err_exit("load ref image");

	if (syn_mode) {

		// Temporary source image
		Image src_temp;	
		init_im(&src_temp);
		if (src_path == NULL) {
			// Copy from the reference
			if (im_copy_data(&ref, &src_temp))
				err_exit("copy reference image\n");
		} else {
			// Load the specified source
			if (read_nii(src_path, &src_temp))
				err_exit("load source image\n");
		}

                // Make a rotation matrix
		double ang_deg = 10.0;
		double ang_rad = ang_deg * UTIL_PI / 180.0;
		MAT_RM_GET(&A, 0, 0, double) = cos(ang_rad);
		MAT_RM_GET(&A, 0, 1, double) = -sin(ang_rad);
		MAT_RM_GET(&A, 1, 0, double) = sin(ang_rad);
		MAT_RM_GET(&A, 1, 1, double) = cos(ang_rad);
		MAT_RM_GET(&A, 2, 2, double) = 1.0;
		if (Affine_set_mat(&A, &aff_syn))
			err_exit("set rotation matrix\n");

                // Rotate the center 
                const double xc = ((double) src_temp.nx - 1.0) / 2.0;
                const double yc = ((double) src_temp.ny - 1.0) / 2.0;
                const double zc = ((double) src_temp.nz - 1.0) / 2.0;
                double xtrans, ytrans, ztrans;
                if (apply_tform_xyz(xc, yc, zc, &xtrans, &ytrans, &ztrans, 
                                        AFFINE, &aff_syn))
                        err_exit("get center translation");

                // Apply the center translation and update the affine transformation
                MAT_RM_GET(&A, 0, 3, double) = xc - xtrans;
                MAT_RM_GET(&A, 1, 3, double) = yc - ytrans;
                MAT_RM_GET(&A, 2, 3, double) = zc - ztrans;
		if (Affine_set_mat(&A, &aff_syn))
			err_exit("set affine matrix\n");

		// Apply an affine transformation to the reference image
		if (im_inv_transform(&src_temp, &src, AFFINE, &aff_syn, INTERP))
			err_exit("apply image transform\n");		

#ifdef VERBOSE
		if (write_Mat_rm(TFORM_IN_PATH, &aff_syn.A))
			err_exit("write input transform matrix to file \n");	
		printf("Input transformation matrix written to %s \n",
			   TFORM_IN_PATH);
		if (print_Mat_rm(&aff_syn.A))
		    err_exit("print input matrix");
#endif

		// Clean up
		im_free(&src_temp);
	} else {
		// Load the source image
		if (read_nii(src_path, &src))
			err_exit("load source image");
	}

	// Zero-pad the images to have the same size
	nx = MAX(src.nx, ref.nx);
	ny = MAX(src.ny, ref.ny);
	nz = MAX(src.nz, ref.nz);
	if (init_im_first_time(&srcp, nx, ny, nz, 1) || 
	    init_im_first_time(&refp, nx, ny, nz, 1) || 
	    im_pad(&src, &srcp) || 
	    im_pad(&ref, &refp))
	    err_exit("pad images \n");

	// Release the original images
	im_free(&src);
	im_free(&ref);

        // Start the benchmark
        printf("Starting the benchmark...\n");
	clock_gettime(CLOCK_MONOTONIC, &reg_start);

	// Extract features from reference image
	if (SIFT3D_detect_keypoints(&detector, &refp, &kp_ref))
		err_exit("detect source keypoints\n");
	if (SIFT3D_extract_descriptors(&extractor, &detector.gpyr, &kp_ref,
		&desc_ref, TRUE))
		err_exit("extract source descriptors\n");


	// Save intermediate data
#ifdef VERBOSE
{
	Image kp_ref_im;
	Mat_rm kp_ref_mat;

	init_im(&kp_ref_im);

	if (write_nii("../debug/src_im.nii", &srcp) ||
		write_nii("../debug/ref_im.nii", &refp))
		err_exit("write input images");
	if (write_pyramid(GPYR_REF_PATH, &detector.gpyr))
		fprintf(stderr, "Failed to write reference GSS pyramid to"
						" path %s \n", GPYR_REF_PATH);
		printf("Reference GSS pyramid written to %s \n", GPYR_REF_PATH);
	if (write_pyramid(DOG_REF_PATH, &detector.dog))
		fprintf(stderr, "Failed to write reference DoG pyramid to"
						" path %s \n", DOG_REF_PATH);
		printf("Reference DoG pyramid written to %s \n", DOG_REF_PATH);
	if (write_Keypoint_store(KP_REF_PATH, &kp_ref))
		err_exit("write keypoints");
	if (init_Mat_rm(&kp_ref_mat, 0, 0, DOUBLE, FALSE) ||
	    Keypoint_store_to_Mat_rm(&kp_ref, &kp_ref_mat) ||
	    draw_points(&kp_ref_mat, refp.nx, refp.ny, refp.nz, 1, &kp_ref_im) ||
	    write_nii(KP_REF_IM_PATH, &kp_ref_im))
	    err_exit("draw reference keypoints");
	printf("%u Reference keypoints written to %s and %s\n",
		   (unsigned int) kp_ref.slab.num, KP_REF_PATH, KP_REF_IM_PATH);

	im_free(&kp_ref_im);
	cleanup_Mat_rm(&kp_ref_mat);
}
#endif

	// Extract source keypoints
	if (SIFT3D_detect_keypoints(&detector, &srcp, &kp_src))
		err_exit("detect source keypoints\n");
	if (SIFT3D_extract_descriptors(&extractor, &detector.gpyr, &kp_src, 
				       &desc_src, TRUE))
		err_exit("extract source descriptors\n");

#ifdef VERBOSE
{
	Image kp_src_im;
	Mat_rm kp_src_mat;

	init_im(&kp_src_im);

	// More intermediate data
	if (write_pyramid(GPYR_SRC_PATH, &detector.gpyr))
		fprintf(stderr, "Failed to write source GSS pyramid to"
				" path %s \n", GPYR_SRC_PATH);
		printf("Source GSS pyramid written to %s \n", GPYR_SRC_PATH);
	if (write_pyramid(DOG_SRC_PATH, &detector.dog))
		fprintf(stderr, "Failed to write source DoG pyramid to"
				" path %s \n", DOG_SRC_PATH);
		printf("Source DoG pyramid written to %s \n", DOG_SRC_PATH);
	if (write_Keypoint_store(KP_SRC_PATH, &kp_src))
		err_exit("write keypoints");
	if (init_Mat_rm(&kp_src_mat, 0, 0, DOUBLE, FALSE) ||
	    Keypoint_store_to_Mat_rm(&kp_src, &kp_src_mat) ||
	    draw_points(&kp_src_mat, srcp.nx, srcp.ny, srcp.nz, 1, &kp_src_im) ||
	    write_nii(KP_SRC_IM_PATH, &kp_src_im))
	    err_exit("draw source keypoints");
	printf("%u Source keypoints written to %s and %s\n",
		   (unsigned int) kp_src.slab.num, KP_SRC_PATH, KP_SRC_IM_PATH);

	im_free(&kp_src_im);
	cleanup_Mat_rm(&kp_src_mat);
}
#endif

	// Match features
	if (SIFT3D_nn_match_fb(&desc_src, &desc_ref, nn_thresh, &matches)) 
		err_exit("match keypoints\n");

	if (SIFT3D_matches_to_Mat_rm(&desc_src, &desc_ref, matches,
				     &match_src, &match_ref))
		err_exit("extract coordinate matrices \n");

#ifdef  VERBOSE
{	
	Image background, overlay;

	init_im(&background);
	init_im(&overlay);

	// Write the matches to a text file
	if (write_Mat_rm(MATCH_SRC_PATH, &match_src) ||	
		write_Mat_rm(MATCH_REF_PATH, &match_ref))
		err_exit("write matches to file\n");

	printf("%u matched features written to %s and %s \n", 
		   (unsigned int) match_src.num_rows, MATCH_SRC_PATH, 
		   MATCH_REF_PATH);
	
	// Draw the matches
	if (draw_matches(&srcp, &refp, &match_src, &match_ref, &background, 
			 &overlay))
		err_exit("draw matches");

	// Save the images
	if (write_nii(OVERLAY_PATH, &overlay) || 
	    write_nii(BG_PATH, &background))
		err_exit("save feature images \n");

	im_free(&background);
	im_free(&overlay);
}
#endif
	
	// Find the transformation 
	if (find_tform_ransac(&ran, &match_src, &match_ref, 3, AFFINE, 
						  (void *) &aff_reg))	
		err_exit("fit transform\n");

	// Transform the source image
	if (im_inv_transform(&srcp, &srcp_reg, AFFINE, &aff_reg, LINEAR))
		err_exit("apply image transform\n");		

	// End the benchmark
	clock_gettime(CLOCK_MONOTONIC, &reg_end);
	diff.tv_sec = reg_end.tv_sec - reg_start.tv_sec;
	diff.tv_nsec = reg_end.tv_nsec - reg_start.tv_nsec;
	elapsed_ms = diff.tv_sec * 1E3 + diff.tv_nsec / 1E6;	
	printf("Registration of source image completed in %ld ms \n", elapsed_ms);

	// Save the result
	if (write_nii(IM_OUT_PATH, &srcp_reg))
		err_exit("write registered image\n");

	if (syn_mode) {

                double precision, recall;
                int i, j, num_false_neg, num_true_pos, num_true_neg;

                const int num_ref = desc_ref.num;
                const int num_src = desc_src.num;
                const int num_pos = match_src.num_rows;
                const int num_neg = num_src - num_pos;

	    // Count the number of true classifications
	    num_true_pos = num_true_neg = 0;
	    for (i = 0; i < num_src; i++) {

		double match_err, min_err, min_err_sq, xst, yst, zst;

		// TODO: Refactor SIFT3D_Descriptor and Keypoint to have Cvec coordinates 

		const SIFT3D_Descriptor *const ds = desc_src.buf + i;
		const double xs = ds->xd;
		const double ys = ds->yd;
		const double zs = ds->zd;

		const int match_idx = matches[i];

		int match;

		// Tranform the source feature by the ground truth
		apply_tform_xyz(xs, ys, zs, &xst, &yst, &zst, AFFINE, &aff_syn);

		// Find the matched feature, if any
		match = match_idx >= 0;

		// Test the match against the ground truth	
		if (match) {
	
		    const SIFT3D_Descriptor *const dr = desc_ref.buf + match_idx;
		    const double mrx = dr->xd;
		    const double mry = dr->yd;
		    const double mrz = dr->zd;

		    // Get the error between the match and the ground truth
		    const double dx = xst - mrx;
		    const double dy = yst - mry;
		    const double dz = zst - mrz;

		    match_err = sqrt(dx * dx + dy * dy + dz * dz);

		    // Test if the match is outside the error threshold
		    if (match_err <= ERR_THRESH)  {
                        // True positive
                        num_true_pos++;
                    }

                    // Nothing left to do for positives
	            continue;
		}
		
		// Find the nearest reference feature location
		min_err_sq = DBL_MAX;
		for (j = 0; j < num_ref; j++) {

		    // Get the error between the reference feature and ground truth
		    const SIFT3D_Descriptor *const dr = desc_ref.buf + j;
		    const double xr = dr->xd;
		    const double yr = dr->yd;
		    const double zr = dr->zd;
		    const double dx = xst - xr;
		    const double dy = yst - yr;
		    const double dz = zst - zr;
		    const double err_sq = dx * dx + dy * dy + dz * dz;

		    // Update the minimum error
		    min_err_sq = MIN(min_err_sq, err_sq);
		}

		min_err = sqrt(min_err_sq);

		if (min_err < ERR_THRESH) {
		    // False negative
		    continue;
                }
   
		// True negative 
		num_true_neg++;
	    }	

	    // Compute the number of false classifications
	    num_false_neg = num_neg - num_true_neg;

	    // Compute the precision and recall
	    precision = (double) num_true_pos / (num_pos + DBL_EPSILON);
	    recall = (double) num_true_pos / (num_true_pos + num_false_neg + 
					      DBL_EPSILON);

	    // Print the statistics
	    printf("\nTotal number of source features: %d \n", num_src);
	    printf("Positives: (%d / %d) true \n", num_true_pos, num_pos);
	    printf("Negatives: (%d / %d) true \n", num_true_neg, num_neg);
	    printf("Precision: %f \nRecall: %f \n\n", precision, recall);
	}

#ifdef VERBOSE
 {
	Image grid, grid_deformed;

	// Write the transformation matrix
	if (write_Mat_rm(TFORM_REG_PATH, &aff_reg.A))
		err_exit("write inferred transform matrix to file\n");	
	printf("Inferred transformation matrix written to %s\n",
			TFORM_REG_PATH);
	if (print_Mat_rm(&aff_reg.A))
	    err_exit("print result matrix \n");

	// Save the deformation grid
	init_im(&grid_deformed);
	if (draw_grid(&grid, nx, ny, nz, 20, 1))
		err_exit("make grid \n");
	if (im_inv_transform(&grid, &grid_deformed, AFFINE, &aff_reg, LINEAR))
		err_exit("deform grid \n");
	if (write_nii(GRID_PATH, &grid_deformed))
		err_exit("write deformation grid to file \n");	
	printf("Deformation grid written to %s \n", GRID_PATH); 
    }
#endif

	return 0;
}
