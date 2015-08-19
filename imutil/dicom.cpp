/* -----------------------------------------------------------------------------
 * dicom.cpp 
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * C-language wrapper for the DCMTK library.
 * -----------------------------------------------------------------------------
 */


/*----------------Include the very picky DCMTK----------------*/
#include "dcmtk/config/osconfig.h"    /* make sure OS specific configuration is included first */

//XXX
#include "dcmtk/config/cfunix.h"

#define INCLUDE_CSTDIO
#define INCLUDE_CSTRING
#include "dcmtk/ofstd/ofstdinc.h"

#ifdef HAVE_GUSI_H
#include <GUSI.h>
#endif

#include "dcmtk/dcmdata/dctk.h"          /* for various dcmdata headers */
#include "dcmtk/dcmdata/cmdlnarg.h"      /* for prepareCmdLineArgs */
#include "dcmtk/dcmdata/dcuid.h"         /* for dcmtk version name */

#include "dcmtk/ofstd/ofconapp.h"        /* for OFConsoleApplication */
#include "dcmtk/ofstd/ofcmdln.h"         /* for OFCommandLine */

#include "dcmtk/oflog/oflog.h"           /* for OFLogger */

#include "dcmtk/dcmimgle/dcmimage.h"     /* for DicomImage */
#include "dcmtk/dcmimage/diregist.h"     /* include to support color images */
#include "dcmtk/dcmdata/dcrledrg.h"      /* for DcmRLEDecoderRegistration */

#ifdef BUILD_DCMSCALE_AS_DCMJSCAL
#include "dcmtk/dcmjpeg/djdecode.h"      /* for dcmjpeg decoders */
#include "dcmtk/dcmjpeg/dipijpeg.h"      /* for dcmimage JPEG plugin */
#endif

#ifdef WITH_ZLIB
#include <zlib.h>          /* for zlibVersion() */
#endif
/*---------------------------------------------------------*/

/* Other includes */
#include <iostream>
#include <stdint.h>
#include "imutil.h"
#include "macros.h"
#include "dicom.h"

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
                std::cerr << tag ": unexpected exception " << std::endl; \
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
        if (dicom.getStatus() != EIS_Normal) {
               std::cerr << "read_dcm_cpp: failed to open image " <<
                        path << " (" << 
                        DicomImage::getString(dicom.getStatus()) << ")" << 
                        std::endl; 
                return SIFT3D_FAILURE;
        }

        // Check for color images
        if (!dicom.isMonochrome()) {
                std::cerr << "read_dcm_cpp: reading of color DICOM " <<
                        "images is not supported at this time" << std::endl;
                return SIFT3D_FAILURE;
        }

        // Set the window 
        dicom.setMinMaxWindow();

        // Get the dimensions
        const int nx = dicom.getWidth();
        const int ny = dicom.getHeight();
        const int nz = dicom.getFrameCount();

        // Resize the output
        im->nx = nx;
        im->ny = ny;
        im->nz = nz;
        im->nc = 1;
        im_default_stride(im);
        if (im_resize(im))
                return SIFT3D_FAILURE;

        // Load the image as a DcmFileFormat 
        DcmFileFormat fileformat;
        OFCondition status = fileformat.loadFile(path);
        if (!status.good()) {
               std::cerr << "read_dcm_cpp: failed to read DICOM file " <<
                        path << " (" << status.text() << ")" << 
                        std::endl; 
                return SIFT3D_FAILURE;
        }

        // Read the pixel spacing
        Float64 pixelSpacing;
        status =  fileformat.getMetaInfo()->findAndGetFloat64(DCM_PixelSpacing, 
                pixelSpacing);
        if (!status.good()) {
                std::cerr << "read_dcm_cpp: failed to get pixel spacing " <<
                        "from file " << path << " (" << status.text() << ")" <<
                        std::endl;
                return SIFT3D_FAILURE;
        }
        const double ux = static_cast<double>(pixelSpacing);
        if (ux <= 0.0) {
                std::cerr << "read_dcm_cpp: file " << path << " has " <<
                        "invalid pixel spacing: " << ux << std::endl;
                return SIFT3D_FAILURE;
        }

        // Read the slice thickness 
        Float64 sliceThickness;
        status = fileformat.getMetaInfo()->findAndGetFloat64(DCM_SliceThickness,
                sliceThickness);
        if (!status.good()) {
                std::cerr << "read_dcm_cpp: failed to get slice thickness " <<
                        "from file " << path << " (" << status.text() << ")" <<
                        std::endl;
        }

        // Convert to double 
        const double uz = static_cast<double>(sliceThickness);
        if (uz  <= 0.0) {
                std::cerr << "read_dcm_cpp: file " << path << " has " <<
                        "invalid slice thickness: " << uz << std::endl;
                return SIFT3D_FAILURE;
        }

        //TODO
        // im->ux = ux;
        // im->uy = ux * dicom.getHeightWidthRatio();
        // im->uz = uz;

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
                                DicomImage::getString(dicom.getStatus()) << 
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
