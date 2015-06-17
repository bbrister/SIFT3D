/* denseSift3d.c 
 * ----------------------------------------------------------------
 * Command-line tool to extract dense SIFT3D features from an image.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "imutil.h"
#include "sift.h"

#define BUF_SIZE (1 << 10)

const char help_msg[] = 
        "Usage: denseSift3d [input.nii] [descriptors%.nii] \n"
        "\n"
        "Extracts a dense gradient histogram image from the input file. The \n"
        "output is a set of 12 images, each representing a channel or \n"
        "histogram bin. The last '%' character in the output filename is \n"
        "replaced by the channel index.\n"
        "\n"
        "Accepted image formats: \n"
        "       -.nii \n"
        "       -.nii.gz \n"
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
void errMsg(const char *msg) {
        fprintf(stderr, "denseSift3d: %s \n", msg);
        fputc('\n', stderr);
        fputs("See \"denseSift3d --help\" for more information. \n", stderr);
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
        if (argc != 3) {
                errMsg("not enough arguments.");
                return 1;
        }
        in_path = argv[1];
        out_path = argv[2];

        /* Initialize data */
        init_im(&im);
        init_im(&desc);
        init_im(&chan);
        if (init_SIFT3D(&sift3d))
                err_exit("initalize sift3d \n");

        /* Read the image */        
        if (read_nii(in_path, &im))
                err_exit("read the image \n");                

        /* Ensure the output file name has a % character */
        if ((marker = strrchr(out_path, '%')) == NULL) {
                errMsg("output filename must contain '%'.");
                return 1;
        }
        marker_pos = marker - out_path;

        /* Get the output file name length */
        len = strlen(out_path) + (int) ceil(log10((double) im.nc)) - 1;
        if (len > BUF_SIZE) {
                fprintf(stderr, "Error: output filename cannot exceed %d "
                                "characters. \n", BUF_SIZE);
                return 1;
        }

        /* Extract the descriptors */
        if (SIFT3D_extract_dense_descriptors(&sift3d, &im, &desc))
                err_exit("extract descriptors \n");
        

        /* Write each channel as a separate image */
        for (c = 0; c < desc.nc; c++) {

                /* Get the channel */
                if (im_channel(&desc, &chan, c))
                        err_exit("extract channel");

                /* Form the output file name */
                out_name[0] = '\0';
                sprintf(chan_str, "%d", c);
                strncat(out_name, out_path, marker_pos);
                strcat(out_name, chan_str);
                strcat(out_name, marker + 1);

                /* Write the channel */
                if (write_nii(out_name, &chan))
                        err_exit(out_path);
        }

        return 0;
}
