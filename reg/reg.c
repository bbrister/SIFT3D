/* -----------------------------------------------------------------------------
 * reg.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This file contains routines for image registration.
 * -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "reg.h"
#include "types.h"
#include "macros.h"
#include "sift.h"
#include "imutil.h"

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
        init_SIFT3D(&reg->sift3d);
        init_im(&reg->src);
        init_im(&reg->ref);
	if (init_Mat_rm(&reg->match_src, 0, 0, DOUBLE, SIFT3D_FALSE) ||
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
        im_free(&reg->src);
        im_free(&reg->ref);
        cleanup_Mat_rm(&reg->match_src);
        cleanup_Mat_rm(&reg->match_ref);
}

/* Set the matching theshold of a Reg_SIFT3D struct. */
int set_nn_thresh_Reg_SIFT3D(Reg_SIFT3D *const reg, const double nn_thresh) {

        if (nn_thresh <= 0) {
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
int set_src_Reg_SIFT3D(Reg_SIFT3D *const reg, Image *const src) {

        SIFT3D *const sift3d = &reg->sift3d;
        Keypoint_store *const kp_src = &reg->kp_src;
        SIFT3D_Descriptor_store *const desc_src = &reg->desc_src;

        // Copy the data
        if (im_copy_data(src, &reg->src))
                return SIFT3D_FAILURE;

        // Detect keypoints
	if (SIFT3D_detect_keypoints(sift3d, &reg->src, kp_src)) {
		fprintf(stderr, "register_SIFT3D: failed to detect source "
                        "keypoints\n");
                return SIFT3D_FAILURE;
        }

        // Extract descriptors
	if (SIFT3D_extract_descriptors(sift3d, kp_src, desc_src)) {
                fprintf(stderr, "register_SIFT3D: failed to extract source "
                                "descriptors \n");
                return SIFT3D_FAILURE;
        }

        return SIFT3D_SUCCESS;
}

/* The same as set_source_Reg_SIFT3D, but sets the reference image. */
int set_ref_Reg_SIFT3D(Reg_SIFT3D *const reg, Image *const ref) {

        SIFT3D *const sift3d = &reg->sift3d;
        Keypoint_store *const kp_ref = &reg->kp_ref;
        SIFT3D_Descriptor_store *const desc_ref = &reg->desc_ref;

        // Copy the data
        if (im_copy_data(ref, &reg->ref))
                return SIFT3D_FAILURE;

        // Detect keypoints
        if (SIFT3D_detect_keypoints(sift3d, &reg->ref, kp_ref)) {
		fprintf(stderr, "register_SIFT3D: failed to detect reference "
                        "keypoints\n");
                return SIFT3D_FAILURE;
        }

        // Extract descriptors
	if (SIFT3D_extract_descriptors(sift3d, kp_ref, desc_ref)) {
		fprintf(stderr, "register_SIFT3D: failed to extract reference "
                        "descriptors\n");
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
	if (find_tform_ransac(ran, match_src, match_ref, 3, tform))
                return SIFT3D_FAILURE;

	return SIFT3D_SUCCESS;
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

