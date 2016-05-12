/* -----------------------------------------------------------------------------
 * denseSift3d.c 
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This file contains a command-line tool to extract dense SIFT3D features from 
 * an image.
 * -----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "immacros.h"
#include "imutil.h"
#include "sift.h"

#define BUF_SIZE (1 << 10)

/* The log tag */
const char tag[] = "denseSift3d";

/* The help message */
const char help_msg[] = 
        "Usage: denseSift3D [input.nii] [descriptors%.nii] \n"
        "\n"
        "Extracts a dense gradient histogram image from the input file. The \n"
        "output is a set of 12 images, each representing a channel or \n"
        "histogram bin. The last '%' character in the output filename is \n"
        "replaced by the channel index.\n"
        "\n"
        "Supported image formats: \n"
	"	.dcm (DICOM) \n"
        "	.nii (nifti-1) \n"
        "	.nii.gz (gzip-compressed nifti-1) \n"
	"	directory containing .dcm files \n"
        "\n"
        "Example: \n"
        "       denseSift3d in.nii.gz out%.nii.gz \n"
        "\n"
        "Upon completion, the output would be the following 12 images: \n"
        "       -out0.nii.gz \n"
        "       -out1.nii.gz \n"
        "            ... \n"
        "       -out11.nii.gz \n"
        "\n";
          
/* Print an error message. */      
void err_msg(const char *msg) {
        SIFT3D_ERR("%s: %s \n"
                "Use \"denseSift3d --help\" for more information. \n", tag, 
                msg);
}

/* Report an unexpected error. */
void err_msgu(const char *msg) {
        err_msg(msg);
        print_bug_msg();
}

int main(int argc, char **argv) {

        char out_name[BUF_SIZE], chan_str[BUF_SIZE];
        Image im, desc, chan;
        SIFT3D sift3d;
        char *in_path, *out_path, *marker;
        size_t len;
        int c, marker_pos;

        /* Parse the GNU standard options */
        switch (parse_gnu(argc, argv)) {
                case SIFT3D_HELP:
                        puts(help_msg);
                        return 0;
                case SIFT3D_VERSION:
                        return 0;
        }

        /* Parse the arguments */
        if (argc < 3) {
                err_msg("Not enough arguments.");
                return 1;
        } else if (argc > 3) {
                err_msg("Too many arguments.");
                return 1;
        }
        in_path = argv[1];
        out_path = argv[2];

        /* Initialize data */
        init_im(&im);
        init_im(&desc);
        init_im(&chan);
        if (init_SIFT3D(&sift3d)) {
                err_msgu("Failed to initialize SIFT3D data.");
                return 1;
        }

        /* Read the image */        
        if (im_read(in_path, &im)) {
                
                char msg[BUF_SIZE];

                snprintf(msg, BUF_SIZE, "Failed to read input image \"%s\".",
			in_path);
                err_msg(msg);
                return 1;
        }

        /* Ensure the output file name has a % character */
        if ((marker = strrchr(out_path, '%')) == NULL) {
                err_msg("output filename must contain '%'.");
                return 1;
        }
        marker_pos = marker - out_path;

        /* Get the output file name length */
        len = strlen(out_path) + (int) ceil(log10((double) im.nc)) - 1;
        if (len > BUF_SIZE) {

                char msg[BUF_SIZE];

                snprintf(msg, BUF_SIZE, "Ouput filename cannot exceed %d "
			"characters.", BUF_SIZE);
                err_msg(msg);
                return 1;
        }

        /* Extract the descriptors */
        if (SIFT3D_extract_dense_descriptors(&sift3d, &im, &desc)) {
                err_msgu("Failed to extract descriptors.");
                return 1;
        }
        

        /* Write each channel as a separate image */
        for (c = 0; c < desc.nc; c++) {

                /* Get the channel */
                if (im_channel(&desc, &chan, c)) {
                        err_msgu("Failed to extract the channel.");
                        return 1;
                }

                /* Form the output file name */
                out_name[0] = '\0';
                snprintf(chan_str, BUF_SIZE, "%d", c);
                strncat(out_name, out_path, marker_pos);
                strcat(out_name, chan_str);
                strcat(out_name, marker + 1);

                /* Write the channel */
                if (im_write(out_name, &chan)) {

                        char msg[BUF_SIZE];

                        snprintf(msg, BUF_SIZE, "Failed to write output image "
				"\"%s\".", out_name);
                        err_msg(msg);
                        return 1;
                }
        }

        return 0;
}
