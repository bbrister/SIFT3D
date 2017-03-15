/* -----------------------------------------------------------------------------
 * regSift3D.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This file contains the CLI to detect and match SIFT3D features between a pair
 * of images, and register one image to another.
 * -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "immacros.h"
#include "imutil.h"
#include "sift.h"
#include "reg.h"

/* Option tags */
#define MATCHES 'a'
#define TRANSFORM 'b'
#define WARPED 'c'
#define CONCAT 'd'
#define KEYS 'e'
#define LINES 'f'
#define NN_THRESH 'g'
#define ERR_THRESH 'h'
#define NUM_ITER 'l'
#define TYPE 'm' 
#define RESAMPLE 'n'

/* Message buffer size */
#define BUF_SIZE 1024

/* Internal parameters */
const interp_type interp = LINEAR; // Interpolation used for the warped image
const tform_type type_default = AFFINE; // Default transformation type
        
/* Print the help message */
static void print_help() {
        printf(
        "Usage: regSift3D [source.nii] [reference.nii] \n"
        "\n"
        "Matches SIFT3D features. \n"
        "\n"
        "Supported input formats: \n"
        " .nii (nifti-1) \n"
        " .nii.gz (gzip-compressed nifti-1) \n"
        "\n"
        "Example: \n"
        " regSift3D --nn_thresh 0.8 --matches matches.csv src.nii ref.nii \n"
        "\n"
        "Output options: \n"
        " --matches [filename] - The feature matches. \n"
        "       Supported file formats: .csv, .csv.gz \n"
        " --transform [filename] - The transformation parameters. \n"
        "       Supported file formats: .csv, .csv.gz \n"
        " --warped [filename] -  The warped source image. \n"
        "       Supported file formats: .dcm, .nii, .nii.gz, directory \n"
        " --concat [filename] - Concatenated images, with source on the left \n"
        "       Supported file formats: .dcm, .nii, .nii.gz, directory \n"
        " --keys [filename] - Keypoints drawn in the concatenated image \n"
        "       Supported file formats: .dcm, .nii, .nii.gz, directory \n"
        " --lines [filename] - Lines drawn between matching keypoints \n"
        "       Supported file formats: .dcm, .nii, .nii.gz, directory \n"
        "At least one output option must be specified. \n"
        "\n"
        "Other options: \n"
        " --nn_thresh [value] - Matching threshold on the nearest neighbor \n"
        "       ratio, in the interval (0, 1]. (default: %.2f) \n"
        " --err_thresh [value] - RANSAC inlier threshold, in the interval \n"
        "       (0, inf). This is a threshold on the squared Euclidean \n"
        "       distance in real-world units. (default: %.1f) \n"
        " --num_iter [value] - Number of RANSAC iterations. (default: %d) \n"
        " --type [value] - Type of transformation to be applied. \n"
        "       Supported arguments: \"affine\" (default: affine) \n"
	" --resample - Internally resample the images to have the same \n"
	"	physical resolution. This is slow. Use it when the images \n"
	"	have very different resolutions, for example registering 5mm \n"
	"	to 1mm slices. \n"
        "\n",
        SIFT3D_nn_thresh_default, SIFT3D_err_thresh_default, 
        SIFT3D_num_iter_default);
        print_opts_SIFT3D();
}

/* Print an error message */
static void err_msg(const char *msg) {
        SIFT3D_ERR("regSift3D: %s \n"
                "Use \"regSift3D --help\" for more information. \n", msg);
}

/* Report an unexpected error. */
static void err_msgu(const char *msg) {
        err_msg(msg);
        print_bug_msg();
}

int main(int argc, char *argv[]) {

        Reg_SIFT3D reg;
        SIFT3D sift3d;
        Ransac ran;
        Image src, ref;
        Mat_rm match_src, match_ref;
        void *tform, *tform_arg;
        char *src_path, *ref_path, *warped_path, *match_path, *tform_path,
                *concat_path, *keys_path, *lines_path;
        tform_type type;
        int num_args, c, have_match, have_tform, resample;

        const struct option longopts[] = {
                {"matches", required_argument, NULL, MATCHES},
                {"transform", required_argument, NULL, TRANSFORM},
                {"warped", required_argument, NULL, WARPED},
                {"concat", required_argument, NULL, CONCAT},
                {"keys", required_argument, NULL, KEYS},
                {"lines", required_argument, NULL, LINES},
                {"nn_thresh", required_argument, NULL, NN_THRESH},
                {"err_thresh", required_argument, NULL, ERR_THRESH},
                {"num_iter", required_argument, NULL, NUM_ITER},
                {"type", required_argument, NULL, TYPE},
		{"resample", no_argument, NULL, RESAMPLE},
                {0, 0, 0, 0}
        };

        const char str_affine[] = "affine";

        // Parse the GNU standard options
        switch (parse_gnu(argc, argv)) {
                case SIFT3D_HELP:
                        print_help();
                        return 0;
                case SIFT3D_VERSION:
                        return 0;
                case SIFT3D_FALSE:
                        break;
                default:
                        err_msgu("Unexpected return from parse_gnu.");
                        return 1;
        }

        // Initialize the data
        init_im(&src);
        init_im(&ref);
        init_Reg_SIFT3D(&reg);
        init_Ransac(&ran);
        if (init_SIFT3D(&sift3d) ||
                init_Mat_rm(&match_src, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE) ||
                init_Mat_rm(&match_ref, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE)) {
                err_msgu("Failed basic initialization.");
                return 1;
        }

        // Initialize parameters to defaults        
        type = type_default;

        // Parse the SIFT3D options
        if ((argc = parse_args_SIFT3D(&sift3d, argc, argv, SIFT3D_FALSE)) < 0)
                return 1;
        if (set_SIFT3D_Reg_SIFT3D(&reg, &sift3d) ||
                set_Ransac_Reg_SIFT3D(&reg, &ran)) {
                err_msgu("Failed to save the SIFT3D or Ransac parameters.");
                return 1;
        }

        // Parse the remaining options 
        opterr = 1;
        have_match = have_tform = resample = SIFT3D_FALSE;
        match_path = tform_path = warped_path = concat_path = keys_path =
                lines_path = NULL;
        while ((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (c) {
                case MATCHES:
                        match_path = optarg;
                        have_match = SIFT3D_TRUE;
                        break;
                case TRANSFORM:
                        tform_path = optarg;
                        have_tform = SIFT3D_TRUE;
                        break;
                case WARPED:
                        warped_path = optarg;
                        have_tform = SIFT3D_TRUE;
                        break;
                case CONCAT:
                        concat_path = optarg;
                        have_match = SIFT3D_TRUE;
                        break;
                case KEYS:
                        keys_path = optarg;
                        have_match = SIFT3D_TRUE;
                        break;
                case LINES:
                        lines_path = optarg;
                        have_match = SIFT3D_TRUE;
                        break;
                case NN_THRESH:
                {
                        const double nn_thresh = atof(optarg);
                        if (set_nn_thresh_Reg_SIFT3D(&reg, nn_thresh)) {
                                err_msg("Invalid value for nn_thresh.");
                                return 1;
                        }
                        break;
                }
                case ERR_THRESH:
                {
                        const double err_thresh = atof(optarg);
                        if (set_err_thresh_Ransac(&ran, err_thresh)) {
                                err_msg("Invalid value for err_thresh.");
                                return 1;
                        }
                        break;
                }
                case NUM_ITER:
                {
                        const int num_iter = atoi(optarg);
                        if (set_num_iter_Ransac(&ran, num_iter)) {
                                err_msg("Invalid value for num_iter.");
                                return 1;
                        }
                        break;
                }
                case TYPE:
                        if (!strcmp(optarg, str_affine)) {
                                type = AFFINE;
                        } else {

                                char msg[BUF_SIZE];

                                snprintf(msg, BUF_SIZE,
                                        "Unrecognized transformation type: %s",
                                        optarg);
                                err_msg(msg);
                                return 1;
                        }
                        break;
                case RESAMPLE:
                        resample = SIFT3D_TRUE;
                        break;
                case '?':
                default:
                        return 1;
                }
        }

        // Ensure that at least one output was specified
        if (!have_match && !have_tform) {
                err_msg("No outputs were specified.");
                return 1;
        }

        // Parse the required arguments
        num_args = argc - optind;
        if (num_args < 2) {
                err_msg("Not enough arguments.");
                return 1;
        } else if (num_args > 2) {
                err_msg("Too many arguments.");
                return 1;
        }
        src_path = argv[optind];
        ref_path = argv[optind + 1];

        // Allocate memory for the transformation
        if ((tform = malloc(tform_type_get_size(type))) == NULL) {
                err_msg("Out of memory.");
                return 1;
        }

        // Initialize the transformation
        if (init_tform(tform, type))
                return 1;

        // Read the images
        if (im_read(src_path, &src)) {

                char msg[BUF_SIZE];

                snprintf(msg, BUF_SIZE, "Failed to read the source image "
			"\"%s\"", src_path);
                err_msg(msg);
                return 1;
        }
        if (im_read(ref_path, &ref)) {

                char msg[BUF_SIZE];

                snprintf(msg, BUF_SIZE, "Failed to read the reference image "
			"\"%s\"", ref_path);
                err_msg(msg);
                return 1;
        }

	// Process the images	
	tform_arg = have_tform ? tform : NULL;
	if (resample) {
		// Optionally register with resampling
		if (register_SIFT3D_resample(&reg, &src, &ref, LINEAR,
			tform_arg)) {
			err_msg("Failed to register the images with "
				 "resampling. \n");
			return 1;
		}
	} else {

		// Set the images
		if (set_src_Reg_SIFT3D(&reg, &src)) {
			err_msg("Failed to set the source image.");
			return 1;
		}
		if (set_ref_Reg_SIFT3D(&reg, &ref)) {
			err_msg("Failed to set the reference image.");
			return 1;
		}

		// Match the features, optionally registering the images 
		if (register_SIFT3D(&reg, tform_arg)) {
			err_msg("Failed to register the images.");
			return 1;
		}
	}

        // Convert the matches to matrices
        if (get_matches_Reg_SIFT3D(&reg, &match_src, &match_ref)) {
                err_msgu("Failed to convert matches to coordinates.");
                return 1;
        }

        // Write the outputs
        if (match_path != NULL) {

                Mat_rm matches;

                // Initialize intermediates
                init_Mat_rm(&matches, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE);

                // Form a combined matrix for both sets of matches
                if (concat_Mat_rm(&match_src, &match_ref, &matches, 1)) {
                        err_msgu("Failed to concatenate the matches.");
                        return 1;
                }

                if (write_Mat_rm(match_path, &matches)) {

                        char msg[BUF_SIZE];

                        snprintf(msg, BUF_SIZE, "Failed to write the matches "
				"\"%s\"", match_path);
                        err_msg(msg);
                        return 1;
                }

                // Clean up
                cleanup_Mat_rm(&matches);
        }
        if (tform_path != NULL && write_tform(tform_path, tform)) {

                char msg[BUF_SIZE];

                snprintf(msg, BUF_SIZE, "Failed to write the transformation "
			"parameters \"%s\"", tform_path);
                err_msg(msg);
                return 1;
        }

        // Optionally warp the source image
        if (warped_path != NULL) {

                Image warped;

                // Initialize intermediates
                init_im(&warped);

                // Set the output to the reference image size
                if (im_copy_dims(&ref, &warped)) {
                        err_msgu("Failed to resize the warped image.");
                        return 1;
                }

                // Warp
                if (im_inv_transform(tform, &src, interp, SIFT3D_FALSE, 
                        &warped)) {
                        err_msgu("Failed to warp the source image.");
                        return 1;
                }

                // Write the warped image
                if (im_write(warped_path, &warped)) {

                        char msg[BUF_SIZE];

                        snprintf(msg, BUF_SIZE, "Failed to write the warped "
				"image \"%s\"", warped_path);
                        err_msg(msg);
                        return 1;
                }

                // Clean up
                im_free(&warped);
        }

        // Optionally draw the matches
        if (concat_path != NULL || keys_path != NULL || lines_path != NULL) {

                Image concat, keys, lines;
                Mat_rm keys_src, keys_ref;

                // Set up the pointers for the images
                Image *const concat_arg = concat_path == NULL ? NULL : &concat;
                Image *const keys_arg = keys_path == NULL ? NULL : &keys;
                Image *const lines_arg = lines_path == NULL ? NULL : &lines;

                // Initialize intermediates
                init_im(&concat);
                init_im(&keys);
                init_im(&lines);
                if (init_Mat_rm(&keys_src, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE) ||
                        init_Mat_rm(&keys_ref, 0, 0, SIFT3D_DOUBLE, 
				SIFT3D_FALSE)) {
                        err_msgu("Failed to initialize keypoint matrices.");
                        return 1;
                }

                // Convert the keypoint coordinates to matrices
                if (SIFT3D_Descriptor_coords_to_Mat_rm(&reg.desc_src, 
                        &keys_src) ||
                    SIFT3D_Descriptor_coords_to_Mat_rm(&reg.desc_ref, 
                        &keys_ref)) {
                        err_msgu("Failed to convert the keypoints to "
                                 "matrices.");
                        return 1;
                }

                // Draw the matches
                if (draw_matches(&src, &ref, &keys_src, &keys_ref, 
                        &match_src, &match_ref, concat_arg, keys_arg,
                        lines_arg)) {
                        err_msgu("Failed to draw the matches.");
                        return 1;
                }

                // Optionally write a concatenated image
                if (concat_path != NULL && im_write(concat_path, &concat)) {

                        char msg[BUF_SIZE];

                        snprintf(msg, BUF_SIZE, "Failed to write the "
				"concatenated image \"%s\"", concat_path);
                        err_msg(msg);
                        return 1;
                }

                // Optionally write the keypoints
                if (keys_path != NULL && im_write(keys_path, &keys)) {

                        char msg[BUF_SIZE];

                        snprintf(msg, BUF_SIZE, "Failed to write the keypoint "
				"image \"%s\"", concat_path);
                        err_msg(msg);
                        return 1;
                }

                // Optionally write the matches
                if (lines_path != NULL && im_write(lines_path, &lines)) {

                        char msg[BUF_SIZE];

                        snprintf(msg, BUF_SIZE, "Failed to write the line "
				"image \"%s\"", concat_path);
                        err_msg(msg);
                        return 1;
                }

                // Clean up
                im_free(&concat);
                im_free(&keys);
                im_free(&lines);
        }

	return 0;
}
