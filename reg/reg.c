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

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "reg.h"
#include "types.h"
#include "macros.h"
#include "sift.h"
#include "imutil.h"

/* Helper routines */
static int subgrad_dist_se(const Image *const ref, const Image *const reg, 
        Image *const grad);
static int subgrad_regular_tich(const Image *const trans, Image *const reg);
static int diffeo_reg_(const Image *const src, const Image *const ref,
                Image *const trans);

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
	if (init_Mat_rm(&reg->match_src, 0, 0, DOUBLE, FALSE) ||
		init_Mat_rm(&reg->match_ref, 0, 0, DOUBLE, FALSE)) {
                fprintf(stderr, "register_SIFT3D: unexpected error \n");
                return FAILURE;
        }

        return SUCCESS;
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
                return FAILURE;
        }
	if (SIFT3D_extract_descriptors(sift3d, (void *const) gpyr,
                                       kp_ref, desc_ref, TRUE)) {
		fprintf(stderr, "register_SIFT3D: failed to extract reference "
                        "descriptors\n");
                return FAILURE;
        }

	// Extract source features 
	if (SIFT3D_detect_keypoints(sift3d, src, kp_src)) {
		fprintf(stderr, "register_SIFT3D: failed to detect source "
                        "keypoints\n");
                return FAILURE;
        }
	if (SIFT3D_extract_descriptors(sift3d, (void *const) gpyr,
                                       kp_src, desc_src, TRUE)) {
                fprintf(stderr, "register_SIFT3D: failed to extract source "
                                "descriptors \n");
                return FAILURE;
        }

	// Match features
	if (SIFT3D_nn_match_fb(desc_src, desc_ref, nn_thresh, matches)) {
		fprintf(stderr, "register_SIFT3D: failed to match "
                                "keypoints \n");
                return FAILURE;
        }

        // Convert matches to coordinate matrices
	if (SIFT3D_matches_to_Mat_rm(desc_src, desc_ref, *matches,
				     match_src, match_ref)) {
		fprintf(stderr, "register_SIFT3D: failed to extract "
                                "coordinate matrices \n");
                return FAILURE;
        }

	// Find the transformation 
	if (find_tform_ransac(ran, match_src, match_ref, 3, type, tform))
		err_exit("fit transform\n");

	return SUCCESS;
}

/* Squared error subgradient. */
static int subgrad_se(const Image *const im, Image *const grad) {

        Cvec grad_vox;
        int x, y, z, c;

        const int x_start = 1;
        const int y_start = 1;
        const int z_start = 1;
        const int x_end = im->nx - 2;
        const int y_end = im->ny - 2;
        const int z_end = im->nz - 2;

        // Resize the gradient
        memcpy(grad->dims, im->dims, IM_NDIMS * sizeof(int));
        grad->nc = IM_NDIMS; 
        im_default_stride(grad);
        if (im_resize(grad))
                return FAILURE;

        // Intialize the gradient to zero
        im_zero(grad);

        IM_LOOP_LIMITED_START(im, x, y, z, x_start, x_end, y_start, y_end,
                z_start, z_end)

                // Iterate over the image channels
                for (c = 0; c < im->nc; c++) {

                        const float vox = IM_GET_VOX(im, x, y, z, c);

                        // Get the image gradient
                        IM_GET_GRAD(im, x, y, z, c, &grad_vox);

                        // Multiply by the image value
                        CVEC_SCALE(&grad_vox, 2 * vox);

                        // Update the gradient
                        IM_GET_VOX(grad, x, y, z, 0) += grad_vox.x;
                        IM_GET_VOX(grad, x, y, z, 1) += grad_vox.y;
                        IM_GET_VOX(grad, x, y, z, 2) += grad_vox.z;
                }
        IM_LOOP_END

        return SUCCESS;
}

/* Squared error distance metric subgradient. */
static int subgrad_dist_se(const Image *const ref, const Image *const reg,
        Image *const grad) {
        
        Image temp1, temp2;

        int x, y, z, c;

        // Allocate the temporary images
        init_im(&temp1);
        init_im(&temp2);
        if (im_copy_dims(ref, &temp1))
                return FAILURE;

        memcpy(temp2.dims, ref->dims, IM_NDIMS * sizeof(int)); 
        temp2.nc = IM_NDIMS;
        im_default_stride(&temp2);
        if (im_resize(&temp2)) {
                im_free(&temp1);
                return FAILURE;
        }

        // Resize grad
        memcpy(grad->dims, ref->dims, IM_NDIMS * sizeof(int));
        grad->nc = IM_NDIMS;
        im_default_stride(grad);
        if (im_resize(grad))
                return FAILURE;

        // Initialize the gradient to zero
        im_zero(grad);

        // Iterate over the image channels
        for (c = 0; c < ref->nc; c++) {

                // Get the difference between reg and ref in this channel
                IM_LOOP_START(&temp1, x, y, z)
                        IM_GET_VOX(&temp1, x, y, z, 0) = 
                                IM_GET_VOX(reg, x, y, z, c) -
                                IM_GET_VOX(ref, x, y, z, c); 
                IM_LOOP_END

                // Get the squared error subgradient for the difference
                if (subgrad_se(&temp1, &temp2))
                        return FAILURE;

                // Add it to the total subgradient
                if (im_acc(&temp2, grad))
                        return FAILURE;
        }

        // Clean up
        im_free(&temp1);
        im_free(&temp2);

        return SUCCESS;
}

/* Tikhonov regularization subgradient. */
static int subgrad_regular_tich(const Image *const trans, Image *const grad) {

        Image temp1, temp2;
        Cvec grad_vox;
        int x, y, z, d;
        float lap;

        const int x_start = 1;
        const int y_start = 1;
        const int z_start = 1;
        const int x_end = trans->nx - 2;
        const int y_end = trans->ny - 2;
        const int z_end = trans->nz - 2;

        // Initialize the temporary images
        init_im(&temp1);
        init_im(&temp2);
        memcpy(temp1.dims, trans->dims, IM_NDIMS * sizeof(int));
        temp1.nc = IM_NDIMS;
        im_default_stride(&temp1);
        if (im_resize(&temp1))
                return FAILURE;
        if (im_copy_dims(&temp1, &temp2)) {
                im_free(&temp1);
                return FAILURE;
        }

        // Resize grad
        memcpy(grad->dims, trans->dims, IM_NDIMS * sizeof(int));
        grad->nc = IM_NDIMS;
        im_default_stride(grad);
        if (im_resize(grad))
                return FAILURE;

        // Initialize the gradient to zero
        im_zero(grad);

        // Iterate over the transform dimensions 
        for (d = 0; d < trans->nc; d++) { 
                // Get the gradient in this transform dimension
                IM_LOOP_LIMITED_START(trans, x, y, z, x_start, x_end, y_start, 
                        y_end, z_start, z_end)

                        IM_GET_GRAD(trans, x, y, z, d, &grad_vox);
                        IM_GET_VOX(&temp1, x, y, z, 0) = grad_vox.x;
                        IM_GET_VOX(&temp1, x, y, z, 1) = grad_vox.y;
                        IM_GET_VOX(&temp1, x, y, z, 2) = grad_vox.z;
                                        
                IM_LOOP_END

                // Get the squared error subgradient
                if (subgrad_se(&temp1, &temp2))
                        return FAILURE;

                // Add it to the total subgradient
                if (im_acc(&temp2, grad))
                        return FAILURE;
        }

        // Clean up
        im_free(&temp1);
        im_free(&temp2);

        return SUCCESS;
}

/* Helper routine to register images with a diffeomorphic transformation. 
 * Starts from an intialized value of trans. */
static int diffeo_reg_(const Image *const src, const Image *const ref,
                Image *const trans) {

        int (*subgrad_dist)(const Image *const, const Image *const, 
                Image *const);
        int (*subgrad_regular)(const Image *const, Image *const);
        Image reg, temp, g_dist, g_regular;
        Image *cur_src, *cur_dst;
        int k;

        // Parameters
        const double alpha = 1E-3;
        const double step = 10;
        const int max_iter = 20;


        // Select the subgradient functions
        subgrad_dist = subgrad_dist_se;
        subgrad_regular = subgrad_regular_tich;

        assert(trans->nc == 3);

        // Initialize the intermediate images
        init_im(&reg);
        init_im(&temp);
        init_im(&g_dist);
        init_im(&g_regular);
        if (im_copy_dims(src, &reg) ||
            im_copy_dims(trans, &temp))
                return FAILURE;

        // Iteratively optimize the transformation
        cur_src = &temp;
        for (k = 0; k < max_iter; k++) {

                int x, y, z, c;

                const int x_start = 1;
                const int y_start = 1;
                const int z_start = 1;
                const int x_end = trans->nx - 2;
                const int y_end = trans->ny - 2;
                const int z_end = trans->nz - 2;
                const int nx = trans->nx;
                const int ny = trans->ny;
                const int nz = trans->nz;

                printf("iteration %d \n", k);

                // Swap the buffers        
                if (cur_src == trans) {
                        cur_src = &temp;
                        cur_dst = trans;
                } else {
                        cur_src = trans;
                        cur_dst = &temp;
                }

                // Transform the source image
                if (im_inv_transform(src, &reg, WARP, cur_src, LINEAR))
                        return FAILURE;

                // Get the distance metric gradient
                if (subgrad_dist(ref, &reg, &g_dist))
                        return FAILURE;

                // Get the regularization gradient
                if (subgrad_regular(cur_src, &g_regular))
                        return FAILURE;

                // Update the transformation
                IM_LOOP_LIMITED_START(trans, x, y, z, x_start, x_end, 
                        y_start, y_end, z_start, z_end)

                        for (c = 0; c < IM_NDIMS; c++) {
                                IM_GET_VOX(cur_dst, x, y, z, c) = 
                                        IM_GET_VOX(cur_src, x, y, z, c) -
                                        step * (IM_GET_VOX(&g_dist, x, y, z, c) + 
                                        alpha * IM_GET_VOX(&g_regular, x, y, z, c));
                        }
#if 0
                        if (x % 50 == 0 && y % 50 == 0 && z % 10 == 0)
                                printf("x: %d y: %d z: %d grad.x: %f \n", x, y, z, 
                                        IM_GET_VOX(cur_dst, x, y, z, 0));
#endif

                        // Clamp to the edges
                        if (IM_GET_VOX(cur_dst, x, y, z, 0) < 0) {
                                IM_GET_VOX(cur_dst, x, y, z, 0) = 0;
                        } else if (IM_GET_VOX(cur_dst, x, y, z, 0) > nx - 1) {
                                IM_GET_VOX(cur_dst, x, y, z, 0) = nx - 1;
                        }
                        if (IM_GET_VOX(cur_dst, x, y, z, 1) < 0) {
                                IM_GET_VOX(cur_dst, x, y, z, 1) = 0;
                        } else if (IM_GET_VOX(cur_dst, x, y, z, 1) > ny - 1) {
                                IM_GET_VOX(cur_dst, x, y, z, 1) = ny - 1;
                        }
                        if (IM_GET_VOX(cur_dst, x, y, z, 2) < 0) {
                                IM_GET_VOX(cur_dst, x, y, z, 2) = 0;
                        } else if (IM_GET_VOX(cur_dst, x, y, z, 2) > nz - 1) {
                                IM_GET_VOX(cur_dst, x, y, z, 2) = nz - 1;
                        }

                IM_LOOP_END

        }

        // Copy the result back to trans, if necessary
        if (cur_dst == &temp && im_copy_data(&temp, trans))
                        return FAILURE;

        // Clean up
        im_free(&temp);
        im_free(&g_dist);
        im_free(&g_regular);
        
        return SUCCESS;
}

/* Register images with a diffeomorphic transformation. */
int diffeo_reg(const Image *const src, const Image *const ref, 
                Image *const trans) {

        Image src_smooth, ref_smooth;
        Gauss_filter gauss;

        double sigma_smooth = 1.0;

        // Verify inputs
        if (memcmp(src->dims, ref->dims, IM_NDIMS * sizeof(int))) {
                fprintf(stderr, "diffeo_reg: images must have the same "
                        "size \n");
                return FAILURE;
        }

        // Initialize the intermediates
        init_im(&src_smooth);
        init_im(&ref_smooth);

        // Intialize the smoothing filter
        if (init_Gauss_filter(&gauss, sigma_smooth, 3))
                return FAILURE;

        // Smooth
        if (apply_Sep_FIR_filter(src, &src_smooth, &gauss.f) ||
            apply_Sep_FIR_filter(ref, &ref_smooth, &gauss.f))
                return FAILURE;

        // Resize the output transformation
	memcpy(trans->dims, src->dims, IM_NDIMS * sizeof(int));	
        trans->nc = IM_NDIMS;
        im_default_stride(trans);
        if (im_resize(trans))
                return FAILURE;

        printf("Source channels: %d reference channels %d \n", src->nc, ref->nc);

        // TODO affine registration

        //TODO: Pyramids

        // Initialize the transformation to zero
        im_zero(trans);

        // Register the images
        if (diffeo_reg_(&src_smooth, &ref_smooth, trans))
                return FAILURE;

        im_free(&src_smooth);
        im_free(&ref_smooth);
        return SUCCESS;
}
