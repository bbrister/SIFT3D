/* -----------------------------------------------------------------------------
 * registerC.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Example of registering two images using the C API.
 */

/* System headers */
#include <stdio.h>

/* SIFT3D headers */
#include "imutil.h"
#include "sift.h"
#include "reg.h"

/* Example file paths */
const char *ref_path = "1.nii.gz";
const char *src_path = "2.nii.gz";
const char *match_path = "1_2_matches.nii.gz";
const char *warped_path = "2_warped.nii.gz";
const char *affine_path = "1_2_affine.csv.gz";

/* This illustrates how to use Reg_SIFT3D within a function, freeing all memory
 * afterwards. */
int demo(void) {

	Image src, ref, warped;
        Reg_SIFT3D reg;
        Affine affine;

        // Initialize the intermediates
        init_im(&src);
        init_im(&ref);
        init_im(&warped);
        if (init_tform(&affine, AFFINE))
                return 1;

        if (init_Reg_SIFT3D(&reg)) {
                cleanup_tform(&affine);
                return 1;
        }

        // Read the images
        if (im_read(src_path, &src) ||
                im_read(ref_path, &ref))
                goto demo_quit;

        // Set the images
        if (set_src_Reg_SIFT3D(&reg, &src) ||
                set_ref_Reg_SIFT3D(&reg, &ref))
                goto demo_quit;

        // Match features and solve for an affine transformation
        if (register_SIFT3D(&reg, &affine))
                goto demo_quit;

        // Write the transformation to a file 
        if (write_tform(affine_path, &affine))
                goto demo_quit;

        // Warp the source image
        if (im_inv_transform(&affine, &src, LINEAR, SIFT3D_TRUE, &warped))
                goto demo_quit;

        // Write the warped image to a file
        if (im_write(warped_path, &warped))
                goto demo_quit;

        // Clean up
        im_free(&src);
        im_free(&ref);
        im_free(&warped);
        cleanup_Reg_SIFT3D(&reg);
        cleanup_tform(&affine);

        return 0;

demo_quit:
        // Clean up and return an error
        im_free(&src);
        im_free(&ref);
        im_free(&warped);
        cleanup_Reg_SIFT3D(&reg);
        cleanup_tform(&affine);

        return 1;
}

int main(void) {

        int ret;

        // Do the demo
        ret = demo();

        // Check for errors
        if (ret != 0) {
                fprintf(stderr, "Fatal demo error, code %d. \n", ret);
                return 1;
        }

        return 0;
}
