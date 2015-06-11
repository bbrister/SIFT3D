/* reg.h
 * ----------------------------------------------------------------
 * Header file for reg.c. 
 * ----------------------------------------------------------------
 */

#include <stdio.h>
#include "reg.h"
#include "types.h"
#include "macros.h"
#include "sift.h"
#include "imutil.h"

/* Initialize a Reg_SIFT3D struct with the default parameters. This must be
 * called before the struct can be used. */
int init_Reg_SIFT3D(Reg_SIFT3D *const reg) {

	reg->matches = NULL;
	init_Keypoint_store(&reg->kp_src);
	init_Keypoint_store(&reg->kp_ref);
	init_SIFT3D_Descriptor_store(&reg->desc_src);
	init_SIFT3D_Descriptor_store(&reg->desc_ref);
	init_Ransac(&reg->ran);
        init_im(&reg->src);
        init_im(&reg->ref);
	if (init_Mat_rm(&reg->match_src, 0, 0, DOUBLE, SIFT3D_FALSE) ||
		init_Mat_rm(&reg->match_ref, 0, 0, DOUBLE, SIFT3D_FALSE)) {
                fprintf(stderr, "register_SIFT3D: unexpected error \n");
                return SIFT3D_FAILURE;
        }

        return SIFT3D_SUCCESS;
}

/* Set the source image. This makes a deep copy of the data, so you are free
 * to modify src after calling this function. */
int set_src_Reg_SIFT3D(Reg_SIFT3D *const reg, Image *const src) {
        return im_copy_data(src, &reg->src);
}

/* The same as set_source_Reg_SIFT3D, but sets the reference image. */
int set_ref_Reg_SIFT3D(Reg_SIFT3D *const reg, Image *const ref) {
        return im_copy_data(ref, &reg->ref);
}

/* Run the registration pipeline. */
int register_SIFT3D(Reg_SIFT3D *const reg, const tform_type type, 
        void *const tform) {

        Ransac *const ran = &reg->ran;
        SIFT3D *const sift3d = &reg->sift3d;
        Pyramid *const gpyr = &sift3d->gpyr;
        Image *const src = &reg->src;
        Image *const ref = &reg->ref;
        Keypoint_store *const kp_src = &reg->kp_src;
        Keypoint_store *const kp_ref = &reg->kp_ref;
        SIFT3D_Descriptor_store *const desc_src = &reg->desc_src;
        SIFT3D_Descriptor_store *const desc_ref = &reg->desc_ref;
        Mat_rm *const match_src = &reg->match_src;
        Mat_rm *const match_ref = &reg->match_ref;
        int **const matches = &reg->matches;
        const double nn_thresh = reg->nn_thresh;

	// Extract reference features 
	if (SIFT3D_detect_keypoints(sift3d, ref, kp_ref)) {
		fprintf(stderr, "register_SIFT3D: failed to detect reference "
                        "keypoints\n");
                return SIFT3D_FAILURE;
        }
	if (SIFT3D_extract_descriptors(sift3d, (void *const) gpyr,
                                       kp_ref, desc_ref, SIFT3D_TRUE)) {
		fprintf(stderr, "register_SIFT3D: failed to extract reference "
                        "descriptors\n");
                return SIFT3D_FAILURE;
        }

	// Extract source features 
	if (SIFT3D_detect_keypoints(sift3d, src, kp_src)) {
		fprintf(stderr, "register_SIFT3D: failed to detect source "
                        "keypoints\n");
                return SIFT3D_FAILURE;
        }
	if (SIFT3D_extract_descriptors(sift3d, (void *const) gpyr,
                                       kp_src, desc_src, SIFT3D_TRUE)) {
                fprintf(stderr, "register_SIFT3D: failed to extract source "
                                "descriptors \n");
                return SIFT3D_FAILURE;
        }

	// Match features
	if (SIFT3D_nn_match_fb(desc_src, desc_ref, nn_thresh, matches)) {
		fprintf(stderr, "register_SIFT3D: failed to match "
                                "keypoints \n");
                return SIFT3D_FAILURE;
        }

        // Convert matches to coordinate matrices
	if (SIFT3D_matches_to_Mat_rm(desc_src, desc_ref, *matches,
				     match_src, match_ref)) {
		fprintf(stderr, "register_SIFT3D: failed to extract "
                                "coordinate matrices \n");
                return SIFT3D_FAILURE;
        }

	// Find the transformation 
	if (find_tform_ransac(ran, match_src, match_ref, 3, type, tform))
		err_exit("fit transform\n");

	return SIFT3D_SUCCESS;
}
