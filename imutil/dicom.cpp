/* -----------------------------------------------------------------------------
 * dicom.cpp 
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * C-language wrapper for the DCMTK library.
 * -----------------------------------------------------------------------------
 */

#include <iostream>
#include <cstdint>

/* Macro to call a C++ function and catch any exceptions it throws,
 * returning SIFT3D_FAILURE when an exception is caught. */
#define CATCH_EXCEPTIONS(tag, fun, ...) \
        try { \
                (fun)( __VA_ARGS__ ); \
        } catch (std::exception &e) { \
\
                std::cerr << e.what() << std::endl; \
\
        } catch (...) { \
\
                std::cerr << TAG ": unexpected exception " << std::endl; \
\
                return SIFT3D_FAILURE; \
        } \

/* Read a DICOM file into an Image struct. */
int read_dcm(const char *path, Image *const im) {

        CATCH_EXCEPTIONS("read_dcm", read_dcm, path, im);

        return SIFT3D_SUCCESS;
}

/* Read all of the DICOM files from a directory into an Image struct. Slices 
 * must be ordered alphanumerically, starting with z = 0. */
int read_dcm_dir(const char *path, Image *const im) {

        CATCH_EXCEPTIONS("read_dcm_dir", read_dcm_dir, path, im);

        return SIFT3D_SUCCESS;
}

/* Write an Image struct into a DICOM file. */
int write_dcm(const char *path, const Image *const im) {

        CATCH_EXCEPTIONS("write_dcm", write_dcm, path, im);

        return SIFT3D_SUCCESS;
}

/* Write an Image struct into a directory of DICOM files. */
int write_dcm_dir(const char *path, const Image *const im) {

        CATCH_EXCEPTIONS("write_dcm_dir", write_dcm_dir, path, im);

        return SIFT3D_SUCCESS;
}

/* Helper function to read a DICOM file using C++ */
static int read_dcm_cpp(const char *path, Image *const im) {

        // Load the image
        DicomImage dicom(path);
        if (image == NULL || dicom.getStatus != EIS_Normal) {
               std::cerr << "read_dcm_dir_cpp: failed to open image " <<
                        path << " (" << 
                        DicomImage::getString(dicom.getStatus()) << ")" << 
                        std::endl; 
                return SIFT3D_FAILURE;
        }

        // Check for color images
        if (!dicom.isMonochrome()) {
                std::cerr << "read_dcm_dir_cpp: reading of color DICOM " <<
                        "images is not supported at this time" << std::endl;
                return SIFT3D_FAILURE;
        }

        // Set the window 
        dicom.setMinMaxWindow();

        // Get the dimensions
        const int nx = getWidth();
        const int ny = getHeight();
        const int nz = dicom.getFrameCount();

        // Resize the output
        im->nx = nx;
        im->ny = ny;
        im->nz = nz;
        im->nc = 1;
        im_default_stride(im);
        if (im_resize(im))
                return SIFT3D_FAILURE;

        // Read each frame
        for (int i = 0; i < nz; i++) { 

                // Get a pointer to the data, rendered as a 32-bit int
                const uint32_t *const frameData = 
                        static_cast<const uint32_t *const >(
                                dicom.getOutputData(32, i));
                if (frameData == NULL) {
                        std::cerr << "read_dcm_dir_cpp: could not get data "
                                << "from image " << path << " frame " << i <<
                                " (" << 
                                DicomImage::getString(image->getStatus()) << 
                                ")" << std::endl; 
                }

                // Copy the frame
                const int x_start = 0;
                const int y_start = 0;
                const int z_start = i;
                const int x_end = nx - 1;
                const int y_end = ny - 1;
                const int z_end = z_start;
                int x, y, z;
                SIFT3D_IM_LOOP_LIMITED_START(im, x, y, z, x_start, x_end, 
                        y_start, y_end, z_start, z_end)

                        SIFT3D_IM_GET_VOX(im, x, y, z, 0) =
                                static_cast<float>(frameData[x + y * nx]);

                SIFT3D_IM_LOOP_END
        }
        
        return SIFT3D_SUCCESS;
}

/* Helper funciton to read a directory of DICOM files using C++ */
static int read_dcm_dir_cpp(const char *path, Image *const im) {

        //TODO: Check if the directory exists, find the dicom images,
        // call read_dcm on them in order, reading to a dummy image,
        // copy the dummy image to each slice

        return SIFT3D_SUCCESS;
}

/* Helper function to write a DICOM file using C++ */
int write_dcm_cpp(const char *path, const Image *const im) {

        //TODO

        return SIFT3D_SUCCESS;
}

/* Helper function to write an image to a directory of DICOM files using C++ */
int write_dcm_dir_cpp(const char *path, const Image *const im) {

        //TODO copy each frame to a dummy image, write it with write_dcm

        return SIFT3D_SUCCESS;
}
