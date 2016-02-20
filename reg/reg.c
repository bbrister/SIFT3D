/* -----------------------------------------------------------------------------
 * reg.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This file contains routines for image registration.
 * -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "reg.h"
#include "types.h"
#include "macros.h"
#include "sift.h"
#include "imutil.h"

/* Internal helper routines */
static void scale_SIFT3D(const double *const factors, 
	Keypoint_store *const kp, SIFT3D_Descriptor_store *const d);

/* Initialize a Reg_SIFT3D struct with the default parameters. This must be
 * called before the struct can be used. */
int init_Reg_SIFT3D(Reg_SIFT3D *const reg) {

	reg->matches = NULL;
        reg->nn_thresh = nn_thresh_default;
	init_Keypoint_store(&reg->kp_src);
	init_Keypoint_store(&reg->kp_ref);
	init_SIFT3D_Descriptor_store(&reg->desc_src);
	init_SIFT3D_Descriptor_store(&reg->desc_ref);
	init_Ransac(&reg->ran);
	if (init_SIFT3D(&reg->sift3d) ||
                init_Mat_rm(&reg->match_src, 0, 0, DOUBLE, SIFT3D_FALSE) ||
		init_Mat_rm(&reg->match_ref, 0, 0, DOUBLE, SIFT3D_FALSE)) {
                fprintf(stderr, "register_SIFT3D: unexpected error \n");
                return SIFT3D_FAILURE;
        }

        return SIFT3D_SUCCESS;
}

/* Free all memory associated with a Reg_SIFT3D struct. reg cannot be reused
 * unless it is reinitialized. */
void cleanup_Reg_SIFT3D(Reg_SIFT3D *const reg) {

        if (reg->matches != NULL)
                free(reg->matches);

        cleanup_Keypoint_store(&reg->kp_src);
        cleanup_Keypoint_store(&reg->kp_ref);
        cleanup_SIFT3D_Descriptor_store(&reg->desc_src);
        cleanup_SIFT3D_Descriptor_store(&reg->desc_ref);
        cleanup_SIFT3D(&reg->sift3d); 
        cleanup_Mat_rm(&reg->match_src);
        cleanup_Mat_rm(&reg->match_ref);
}

/* Set the matching theshold of a Reg_SIFT3D struct. */
int set_nn_thresh_Reg_SIFT3D(Reg_SIFT3D *const reg, const double nn_thresh) {

        if (nn_thresh <= 0 || nn_thresh > 1) {
                fprintf(stderr, "set_nn_thresh_Reg_SIFT3D: invalid threshold: "
                        "%f \n", nn_thresh);
                return SIFT3D_FAILURE;
        }

        reg->nn_thresh = nn_thresh;
        return SIFT3D_SUCCESS;
}

/* Set the SIFT3D parameters of the Reg_SIFT3D struct. Makes a deep copy of
 * sift3d, so you are free to modify it after calling this function. */
void set_SIFT3D_Reg_SIFT3D(Reg_SIFT3D *const reg, const SIFT3D *const sift3d) {
        copy_SIFT3D(sift3d, &reg->sift3d);
}

/* Set the source image. This makes a deep copy of the data, so you are free
 * to modify src after calling this function. */
int set_src_Reg_SIFT3D(Reg_SIFT3D *const reg, const Image *const src) {

        SIFT3D *const sift3d = &reg->sift3d;
        Keypoint_store *const kp_src = &reg->kp_src;
        SIFT3D_Descriptor_store *const desc_src = &reg->desc_src;

        // Detect keypoints
	if (SIFT3D_detect_keypoints(sift3d, src, kp_src)) {
		fprintf(stderr, "set_src_Reg_SIFT3D: failed to detect source "
                        "keypoints\n");
                return SIFT3D_FAILURE;
        }

        // Extract descriptors
	if (SIFT3D_extract_descriptors(sift3d, kp_src, desc_src)) {
                fprintf(stderr, "set_ref_Reg_SIFT3D: failed to extract source "
                                "descriptors \n");
                return SIFT3D_FAILURE;
        }

        return SIFT3D_SUCCESS;
}

/* The same as set_source_Reg_SIFT3D, but sets the reference image. */
int set_ref_Reg_SIFT3D(Reg_SIFT3D *const reg, const Image *const ref) {

        SIFT3D *const sift3d = &reg->sift3d;
        Keypoint_store *const kp_ref = &reg->kp_ref;
        SIFT3D_Descriptor_store *const desc_ref = &reg->desc_ref;

        // Detect keypoints
        if (SIFT3D_detect_keypoints(sift3d, ref, kp_ref)) {
		fprintf(stderr, "set_ref_Reg_SIFT3D: failed to detect "
			"reference keypoints\n");
                return SIFT3D_FAILURE;
        }

        // Extract descriptors
	if (SIFT3D_extract_descriptors(sift3d, kp_ref, desc_ref)) {
		fprintf(stderr, "set_ref_Reg_SIFT3D: failed to extract "
			"reference descriptors\n");
                return SIFT3D_FAILURE;
        }

        return SIFT3D_SUCCESS;
}

/* Run the registration procedure. 
 *
 * Parameters: 
 *   reg: The struct holding registration state.
 *   tform: The output transformation. If NULL, this function only performs
 *     feature matching.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int register_SIFT3D(Reg_SIFT3D *const reg, void *const tform) {

        Ransac *const ran = &reg->ran;
        Mat_rm *const match_src = &reg->match_src;
        Mat_rm *const match_ref = &reg->match_ref;
        int **const matches = &reg->matches;
        const double nn_thresh = reg->nn_thresh;
        SIFT3D_Descriptor_store *const desc_src = &reg->desc_src;
        SIFT3D_Descriptor_store *const desc_ref = &reg->desc_ref;

	// Verify inputs
	if (desc_src->num <= 0) {
		fprintf(stderr, "register_SIFT3D: no source image descriptors "
			"are available \n");
		return SIFT3D_FAILURE;
	}
	if (desc_ref->num <= 0) {
		fprintf(stderr, "register_SIFT3D: no reference image "
			"descriptors are available \n");
		return SIFT3D_FAILURE;
	}

	// Match features
	if (SIFT3D_nn_match_fb(desc_src, desc_ref, nn_thresh, matches)) {
		fprintf(stderr, "register_SIFT3D: failed to match "
                                "descriptors \n");
                return SIFT3D_FAILURE;
        }

        // Convert matches to coordinate matrices
	if (SIFT3D_matches_to_Mat_rm(desc_src, desc_ref, *matches,
				     match_src, match_ref)) {
		fprintf(stderr, "register_SIFT3D: failed to extract "
                                "coordinate matrices \n");
                return SIFT3D_FAILURE;
        }

        // Quit if no tform was provided
        if (tform == NULL)
                return SIFT3D_SUCCESS;

	// Find the transformation 
	if (find_tform_ransac(ran, match_src, match_ref, tform))
                return SIFT3D_FAILURE;

	return SIFT3D_SUCCESS;
}

/* Helper function to scale the keypoints and descriptors by the given 
 * factors */
static void scale_SIFT3D(const double *const factors, 
	Keypoint_store *const kp, SIFT3D_Descriptor_store *const d) {

	double det, scale_factor;
	int i, j, k;

	// Compute the determinant of the scaling transformation
	det = 1.0;
	for (i = 0; i < IM_NDIMS; i++) {
		det *= factors[i];	
	}

	// Compute the scale parameter factor from the determinant
	scale_factor = pow(det, -1.0 / (double) IM_NDIMS);
	
	// Scale the keypoints
	for (k = 0; k < kp->slab.num; k++) {

		Keypoint *const key = kp->buf + k;
		Mat_rm *const R = &key->R;

		// Scale the coordinates
		key->xd *= factors[0];
		key->yd *= factors[1];
		key->zd *= factors[2];

		// Adjust the scale parameter
		key->sd *= scale_factor;

		// Adjust the orientation matrix
		SIFT3D_MAT_RM_LOOP_START(R, i, j)
			SIFT3D_MAT_RM_GET(R, i, j, float) *= 
				(float) (factors[j] / det);
		SIFT3D_MAT_RM_LOOP_END
	}

	// Scale the descriptors
	for (i = 0; i < d->num; i++) {

		SIFT3D_Descriptor *const desc = d->buf + i;

		// Scale the coordinates
		desc->xd *= factors[0];
		desc->yd *= factors[1];
		desc->zd *= factors[2];

		// Adjust the scale parameter
		desc->sd *= scale_factor;
	}
}

/* Like register_SIFT3D, but resamples the input images to have the same
 * physical resolution before extracting features. Use this when registering 
 * images with very different resolutions. The results are converted to the
 * original resolution.
 *
 * Parameters:
 *   reg: See register_SIFT3D.
 *   src: The source, or moving image.
 *   ref: The reference, or fixed image.
 *   interp: The type of interpolation to use.
 *   tform: See register_SIFT3D. 
 *
 * Note that some fields of the returned keypoints, such as scale and
 * orientation, will not make sense in the new coordinate system.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int register_SIFT3D_resample(Reg_SIFT3D *const reg, const Image *const src,
	const Image *const ref, const interp_type interp, void *const tform) {

	double units_min[IM_NDIMS], factors_src[IM_NDIMS], 
		factors_ref[IM_NDIMS];
	Image src_interp, ref_interp;
	int i;

	// Check for the trivial case, when src and dst have the same units
	if (!memcmp(SIFT3D_IM_GET_UNITS(src), SIFT3D_IM_GET_UNITS(ref), 
		IM_NDIMS * sizeof(double))) {
		return set_src_Reg_SIFT3D(reg, src) ||
			set_ref_Reg_SIFT3D(reg, ref) ||
			register_SIFT3D(reg, tform) ? 
			SIFT3D_FAILURE : SIFT3D_SUCCESS;
	}

	// Initalize intermediates
	init_im(&src_interp);
	init_im(&ref_interp);

	// Compute the new units and scaling factors
	for (i = 0; i < IM_NDIMS; i++) {
		const double unit_src = SIFT3D_IM_GET_UNITS(src)[i];
		const double unit_ref = SIFT3D_IM_GET_UNITS(ref)[i];

		// Compute the minimum units between the two images
		units_min[i] = SIFT3D_MIN(unit_src, unit_ref);

		// Compute the scaling factors between the interpolated
		// images and the originals
		factors_src[i] = units_min[i] / unit_src;
		factors_ref[i] = units_min[i] / unit_ref;
	}

	// Resample the images
	if (im_resample(src, units_min, interp, &src_interp) ||
		im_resample(ref, units_min, interp, &ref_interp))
		goto register_interp_quit;

	// Extract features from the interpolated images
	if (set_src_Reg_SIFT3D(reg, &src_interp) ||
		set_ref_Reg_SIFT3D(reg, &ref_interp))
		goto register_interp_quit;

	// Convert the keypoints and descriptors to the original units
	scale_SIFT3D(factors_src, &reg->kp_src, &reg->desc_src);
	scale_SIFT3D(factors_ref, &reg->kp_ref, &reg->desc_ref);

	// Register the images
	if (register_SIFT3D(reg, tform))
		goto register_interp_quit;

	// Clean up
	im_free(&src_interp);	
	im_free(&ref_interp);	

	return SIFT3D_SUCCESS;

register_interp_quit:
	im_free(&src_interp);	
	im_free(&ref_interp);	
	return SIFT3D_FAILURE;
}

/* Write the coordinates of matching keypoints to the matrices match_src
 * and match_ref. This function uses the keypoints and the matches from
 * the last call to register_SIFT3D() on this Reg_SIFT3D struct. */
int get_matches_Reg_SIFT3D(const Reg_SIFT3D *const reg, Mat_rm *const match_src,
        Mat_rm *const match_ref) {

        // Check if we have any matches
        if (reg->matches == NULL)
                return SIFT3D_FAILURE;

        // Copy the matches
        return copy_Mat_rm(&reg->match_src, match_src) ||
                copy_Mat_rm(&reg->match_ref, match_ref);
}

