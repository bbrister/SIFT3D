/* -----------------------------------------------------------------------------
 * ioC.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Short example of file I/O, converting a NIFTI file to DICOM.
 */

/* System headers */
#include <stdio.h>

/* SIFT3D headers */
#include "imutil.h"

/* Example file paths */
const char in_path[] = "1.nii.gz";
const char out_path[] = "1dicom";

/* This illustrates how to use images within a function, and free all memory
 * afterwards. */
int demo(void) {

	Image im;

        // You must call this before any Image struct can be used
        init_im(&im);

        // Read the NIFTI image as input
        if (im_read(in_path, &im))
                goto demo_quit;

        // Write it to a DICOM series
        if (im_write(out_path, &im))
                goto demo_quit;

        // Clean up
        im_free(&im);

        return 0;

demo_quit:
        // Clean up and return an error
        im_free(&im);
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
