/* matchSift3D.c
 * ----------------------------------------------------------------
 * This file contains the CLI to detect and match SIFT3D features. */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "types.h"
#include "imutil.h"
#include "sift.h"
#include "reg.h"

/* Option tags */
#define MATCHES 'a'
#define TRANSFORM 'b'
#define WARPED 'c'
#define THRESH 'd'
#define TYPE 'e' 

/* The help message */
const char help_msg[] = 
        "Usage: regSift3D [source.nii] [reference.nii] \n"
        "\n"
        "Matches SIFT3D features. \n"
        "\n"
        "Supported input formats: \n"
        " .nii (nifti-1) \n"
        " .nii.gz (gzip-compressed nifti-1) \n"
        "\n"
        "Example: \n"
        " regSift3D --thresh 0.8 --matches matches.csv im1.nii im2.nii \n"
        "\n"
        "Output options: \n"
        " --matches [filename] - Writes the feature matches. \n"
        "       Supported file formats: .csv, .csv.gz \n"
        " --transform [filename] - Writes the transformation parameters. \n"
        "       Supported file formats: .csv, .csv.gz \n"
        " --warped [filename] -  Writes the warped image. \n"
        "       Supported file formats: .nii, .nii.gz \n"
        "At least one output option must be specified. \n"
        "\n"
        "Other options: \n"
        " --thresh [value] - Matching threshold on the nearest neighbor \n"
        "       ratio, in the interval (0, 1]. (default: %f) \n"
        " --type [value] - Type of transformation to be applied. \n"
        "       Supported arguments: \"affine\" (default: affine) \n"
        "\n";

/* External parameters */
extern const double nn_thresh_default; // Default matching threshold

/* Internal parameters */
const interp_type interp = LINEAR; // Interpolation used for the warped image
const tform_type type_default = AFFINE; // Default transformation type
        
/* Print the help message */
void print_help() {
        printf(help_msg, nn_thresh_default);
        print_opts_SIFT3D();
}

/* Print an error message */
void err_msg(const char *msg) {
        fprintf(stderr, "regSift3D: %s \n"
                "Use \"regSift3D --help\" for more information. \n", msg);
}

/* Report an unexpected error. */
void err_msgu(const char *msg) {
        err_msg(msg);
        print_bug_msg();
}

int main(int argc, char *argv[]) {

        SIFT3D sift3d;
        Reg_SIFT3D reg;
        Image src, ref;
        void *tform;
        char *src_path, *ref_path, *warped_path, *match_path, *tform_path;
        tform_type type;
        double thresh;
        int num_args, c;

        const struct option longopts[] = {
                {"matches", required_argument, NULL, MATCHES},
                {"transform", required_argument, NULL, TRANSFORM},
                {"warped", required_argument, NULL, WARPED},
                {"thresh", required_argument, NULL, THRESH},
                {"type", required_argument, NULL, TYPE},
                {0, 0, 0, 0}
        };

        const char str_affine[] = "affine";

        // Parse the GNU standard options
        switch (parse_gnu(argc, argv)) {
                case SIFT3D_HELP:
                        puts(help_msg);
                        print_opts_SIFT3D();
                        return 0;
                case SIFT3D_VERSION:
                        return 0;
        }

        // Initialize the data
        init_im(&src);
        init_im(&ref);
        init_Reg_SIFT3D(&reg);
        init_SIFT3D(&sift3d);

        // Initialize parameters to defaults        
        tform = NULL;
        warped_path = match_path = tform_path = NULL;
        type = type_default;

        // Parse the SIFT3D options
        parse_args_SIFT3D(&sift3d, argc, argv, &optind, 0);
        set_SIFT3D_Reg_SIFT3D(&reg, &sift3d);

        // Parse the remaining options 
        opterr = 1;
        while ((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (c) {
                        case MATCHES:
                                match_path = optarg;
                                break;
                        case TRANSFORM:
                                tform_path = optarg;
                                break;
                        case WARPED:
                                warped_path = optarg;
                                break;
                        case THRESH:
                                thresh = atof(optarg);
                                if (set_nn_thresh_Reg_SIFT3D(&reg, thresh)) {
                                        err_msg("Invalid value for thresh.");
                                        return 1;
                                }
                                break;
                        case TYPE:
                                if (!strcmp(optarg, str_affine)) {
                                        type = AFFINE;
                                } else {

                                        char msg[1024];

                                        sprintf(msg,    
                                                "Unrecognized transformation "
                                                "type: %s\n", optarg);
                                        err_msg(msg);
                                        return 1;
                                }
                                break;
                        case '?':
                        default:
                                return 1;
                }
        }

        // Ensure that at least one output was specified
        if (match_path == NULL && tform_path == NULL && warped_path == NULL) {
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
        if (read_nii(src_path, &src)) {

                char msg[1024];

                sprintf(msg, "Failed to read the source image \"%s\"", 
                        src_path);
                err_msg(msg);
                return 1;
        }
        if (read_nii(ref_path, &ref)) {

                char msg[1024];

                sprintf(msg, "Failed to read the reference image \"%s\"",
                        ref_path);
                err_msg(msg);
                return 1;
        }

        // Set the images
        if (set_src_Reg_SIFT3D(&reg, &src)) {
                err_msgu("Failed to set the source image.");
                return 1;
        }
        if (set_ref_Reg_SIFT3D(&reg, &ref)) {
                err_msgu("Failed to set the reference image.");
                return 1;
        }

        // Register
        if (register_SIFT3D(&reg, tform)) {
                err_msgu("Failed to register the images.");
                return 1;
        }

        // Write the outputs
        if (match_path != NULL) {

                Mat_rm matches, match_src, match_ref;
                int i, j;

                SIFT3D_Descriptor_store *const desc_src = &reg.desc_src;
                SIFT3D_Descriptor_store *const desc_ref = &reg.desc_ref;

                // Initialize intermediates
                init_Mat_rm(&matches, 0, 0, DOUBLE, SIFT3D_FALSE);
                init_Mat_rm(&match_src, 0, 0, DOUBLE, SIFT3D_FALSE);
                init_Mat_rm(&match_ref, 0, 0, DOUBLE, SIFT3D_FALSE);

                // Convert the matches to matrices
	        if (SIFT3D_matches_to_Mat_rm(desc_src, desc_ref, reg.matches,
		        &match_src, &match_ref)) {
                        err_msgu("Failed to convert matches to coordinates.");
                        return 1;
                }

                // Form a combined matrix for both sets of matches
                if (concat_h_Mat_rm(&match_src, &match_ref, &matches)) {
                        err_msgu("Failed to concatenate the matches.");
                        return 1;
                }

                if (write_Mat_rm(match_path, &matches)) {

                        char msg[1024];

                        sprintf(msg, "Failed to write the matches \"%s\"", 
                                match_path);
                        err_msg(msg);
                        return 1;
                }

                // Clean up
                cleanup_Mat_rm(&matches);
                cleanup_Mat_rm(&match_src);
                cleanup_Mat_rm(&match_ref);
        }
        if (tform_path != NULL && write_tform(tform_path, tform)) {

                char msg[1024];

                sprintf(msg, "Failed to write the transformation parameters "
                        "\"%s\"", tform_path);
                err_msg(msg);
                return 1;
        }

        // Optionally warp the source image
        if (warped_path != NULL) {

                Image warped;

                // Initialize intermediates
                init_im(&warped);

                // Warp
                if (im_inv_transform(tform, &src, &warped, interp)) {
                        err_msgu("Failed to warp the source image.");
                        return 1;
                }

                // Write the warped image
                if (write_nii(warped_path, &warped)) {

                        char msg[1024];

                        sprintf(msg, "Failed to write the warped image \"%s\"",
                                warped_path);
                        err_msg(msg);
                        return 1;
                }

                // Clean up
                im_free(&warped);
        }

        //TODO: add drawOverlay, drawMatches functions to reg.c, allow them as
        // outputs

	return 0;
}
