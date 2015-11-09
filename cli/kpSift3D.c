/* -----------------------------------------------------------------------------
 * kpSift3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This file contains the CLI to extract SIFT3D keypoints and descriptors from
 * a single image. 
 * -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <getopt.h>
#include "imutil.h"
#include "sift.h"

/* Options */
#define KEYS 'a'
#define DESC 'b'
#define DRAW 'c'

/* Help message */
const char help_msg[] = 
        "Usage: kpSift3D [image.nii] \n"
        "\n"
        "Detects SIFT3D keypoints and extracts their descriptors from an "
        "image.\n" 
        "\n"
        "Example: \n"
        " kpSift3D --keys keys.csv --desc desc.csv image.nii \n"
        "\n"
        "Output options: \n"
        " --keys [filename] \n"
        "       Specifies the output file name for the keypoints. \n"
        "       Supported file formats: .csv, .csv.gz \n"
        " --desc [filename] \n"
        "       Specifies the output file name for the descriptors. \n"
        "       Supported file formats: .csv, .csv.gz \n"
        " --draw [filename] \n"
        "       Draws the keypoints in image space. \n"
        "       Supported file formats: .nii, .nii.gz \n"
        "At least one of the output options must be specified. \n"
        "\n";

/* Print an error message */
void err_msg(const char *msg) {
        fprintf(stderr, "kpSift3D: %s \n"
                "Use \"kpSift3D --help\" for more information. \n", msg);
}

/* Report an unexpected error. */
void err_msgu(const char *msg) {
        err_msg(msg);
        print_bug_msg();
}

/* CLI for 3D SIFT */
int main(int argc, char *argv[]) {

	Image im;
	SIFT3D sift3d;
	Keypoint_store kp;
	SIFT3D_Descriptor_store desc;
	char *im_path, *keys_path, *desc_path, *draw_path;
        int c, num_args;

        const struct option longopts[] = {
                {"keys", required_argument, NULL, KEYS},
                {"desc", required_argument, NULL, DESC},
                {"draw", required_argument, NULL, DRAW},
                {0, 0, 0, 0}
        };

        // Parse the GNU standard options
        switch (parse_gnu(argc, argv)) {
                case SIFT3D_HELP:
                        puts(help_msg);
                        print_opts_SIFT3D();
                        return 0;
                case SIFT3D_VERSION:
                        return 0;
                case SIFT3D_FALSE:
                        break;
                default:
                        err_msgu("Unexpected return from parse_gnu \n");
                        return 1;
        }

	// Initialize the SIFT data 
	if (init_SIFT3D(&sift3d)) {
		err_msgu("Failed to initialize SIFT data.");
                return 1;
        }

        // Parse the SIFT3D options and increment the argument list
        if ((argc = parse_args_SIFT3D(&sift3d, argc, argv, SIFT3D_FALSE)) < 0)
                return 1;

        // Parse the kpSift3d options
        opterr = 1;
        keys_path = desc_path = draw_path = NULL;
        while ((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (c) {
                        case KEYS:
                                keys_path = optarg;
                                break;
                        case DESC:
                                desc_path = optarg;
                                break;
                        case DRAW:
                                draw_path = optarg;
                                break;
                        case '?':
                        default:
                                return 1;
                }
        }

        // Ensure we have at least one output
        if (keys_path == NULL && desc_path == NULL && draw_path == NULL) {
                err_msg("No outputs specified.");
                return 1;
        }

        // Parse the required arguments
        num_args = argc - optind;
        if (num_args < 1) {
                err_msg("Not enough arguments.");
                return 1;
        } else if (num_args > 1) {
                err_msg("Too many arguments.");
                return 1;
        }
        im_path = argv[optind];

	// Initialize data 
	init_Keypoint_store(&kp); 
	init_SIFT3D_Descriptor_store(&desc); 
	init_im(&im);

	// Read the image
	if (im_read(im_path, &im)) {
		err_msg("Could not read image.");
                return 1;
        }

	// Extract keypoints
	if (SIFT3D_detect_keypoints(&sift3d, &im, &kp)) {
		err_msgu("Failed to detect keypoints.");
                return 1;
        }

        // Optionally write the keypoints 
        if (keys_path != NULL && write_Keypoint_store(keys_path, &kp)) {

                char msg[1024];

                sprintf(msg, "Failed to write the keypoints to \"%s\"", 
                        keys_path);
                err_msg(msg);
                return 1;
        }

        // Optionally extract descriptors
        if (desc_path != NULL) {

                // Extract descriptors
	        if (SIFT3D_extract_descriptors(&sift3d, &kp,&desc)) {
                        err_msgu("Failed to extract descriptors.");
                        return 1;
                }

                // Write the descriptors
                if (write_SIFT3D_Descriptor_store(desc_path, &desc)) {

                        char msg[1024];

                        sprintf(msg, "Failed to write the descriptors to "
                                "\"%s\"", desc_path);
                        err_msg(msg);
                        return 1;
                }
        }

        // Optionally draw the keypoints
        if (draw_path != NULL) {

                Image draw;
                Mat_rm keys;

                // Initialize intermediates
                init_im(&draw);
                if (init_Mat_rm(&keys, 0, 0, DOUBLE, SIFT3D_FALSE))
                        err_msgu("Failed to initialize keys matrix");

                // Convert to matrices
                if (Keypoint_store_to_Mat_rm(&kp, &keys)) {
                        err_msgu("Failed to convert the keypoints to "
                                 "a matrix.");
                        return 1;
                }

                // Draw the points
                if (draw_points(&keys, im.dims, 1, &draw)) {
                        err_msgu("Failed to draw the points.");
                        return 1;
                }

                // Write the output
                if (im_write(draw_path, &draw)) {
                        
                        char msg[1024];

                        sprintf(msg, "Failed to draw the keypoints to \"%s\"",
                                draw_path);
                        err_msg(msg);
                        return 1;
                }

                // Clean up
                im_free(&draw);
        }

	return 0;
}
