/* kpSift3D.c
 * ----------------------------------------------------------------
 * This file contains the CLI to run SIFT3D on a volumetric image. It 
 * calls the routines in sift.c.
 */

#include <stdio.h>
#include <getopt.h>
#include "sift.h"

/* Options */
#define KEYS 'a'
#define DESC 'b'

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
        " --desc [filename] \n"
        "       Specifies the output file name for the descriptors. \n"
        "At least one of the output options must be specified. \n"
        "\n"
        "Supported input formats: \n"
        " .nii (nifti-1) \n"
        " .nii.gz (gzip-compressed nifti-1) \n"
        "\n"
        "Supported output formats: \n"
        " .csv (comma-separated value) \n"
        " .csv.gz (gzip-compressed comma-separated value) \n"
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
	char *im_path, *keys_path, *desc_path;
        int c, num_args;

        const struct option longopts[] = {
                {"keys", required_argument, NULL, KEYS},
                {"desc", required_argument, NULL, DESC},
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
        }

	// Initialize the SIFT data 
	if (init_SIFT3D(&sift3d)) {
		err_msgu("Failed to initialize SIFT data.");
                return 1;
        }

        // Parse the SIFT3D options and increment the argument list
        parse_args_SIFT3D(&sift3d, argc, argv, &optind, 0);

        // Parse the kpSift3d options
        opterr = 1;
        keys_path = desc_path = NULL;
        while ((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (c) {
                        case KEYS:
                                keys_path = optarg;
                                //TODO: check the file extension
                                break;
                        case DESC:
                                desc_path = optarg;
                                //TODO: check the file extension
                                break;
                        case '?':
                        default:
                                return 1;
                }
        }

        // Ensure we have at least one output
        if (keys_path == NULL && desc_path == NULL) {
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
        im_path = argv[optind]; //TODO: check the file extension

	// Initialize data 
	init_Keypoint_store(&kp); 
	init_SIFT3D_Descriptor_store(&desc); 
	init_im(&im);

	// Read the image
	if (read_nii(im_path, &im)) {
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

        // Extract descriptors
        if (desc_path == NULL) {

                return 0;

	} else if (SIFT3D_extract_descriptors(&sift3d, &sift3d.gpyr, &kp,
		&desc, SIFT3D_TRUE)) {
		err_msgu("Failed to extract descriptors.");
                return 1;
        }

        // Write the descriptors
        if (write_SIFT3D_Descriptor_store(desc_path, &desc)) {

                char msg[1024];

                sprintf(msg, "Failed to write the descriptors to \"%s\"",
                        desc_path);
                err_msg(msg);
                return 1;
        }

	return 0;
}
