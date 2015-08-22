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

#include "dcmtk/dcmjpeg/djrplol.h"  /* for DJ_RPLossless */

#include "dcmtk/dcmsr/dsrdoc.h" /* DSR report handling */

#ifdef WITH_ZLIB
#include <zlib.h>          /* for zlibVersion() */
#endif
/*---------------------------------------------------------*/

/* Other includes */
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>
#include <stdint.h>
#include <dirent.h>
#include "imutil.h"
#include "macros.h"
#include "dicom.h"

/* Macro to call a C++ function and catch any exceptions it throws,
 * returning SIFT3D_FAILURE when an exception is caught. The return value is
 * stored in ret. */
#define CATCH_EXCEPTIONS(ret, tag, fun, ...) \
        try { \
                ret = (fun)( __VA_ARGS__ ); \
        } catch (std::exception &e) { \
                std::cerr << tag << ": " << e.what() << std::endl; \
                ret = SIFT3D_FAILURE; \
        } catch (...) { \
                std::cerr << tag << ": unexpected exception " << std::endl; \
                ret = SIFT3D_FAILURE; \
        } \

/* File separator */
#ifdef _WINDOWS
const char sep = '\\';
#else
const char sep = '/';
#endif

/* Helper declarations */
static int read_dcm_cpp(const char *path, Image *const im);
static int read_dcm_dir_cpp(const char *path, Image *const im);
static int write_dcm_cpp(const char *path, const Image *const im);
static int write_dcm_dir_cpp(const char *path, const Image *const im);

/* Helper class to store DICOM data. */
class Dicom {
private:
        std::string filename; // DICOM file name
        long long int patient; // Patient ID
        long long int study; // Study ID 
        long long int series; // Series number
        long long int instance; // Identifying instance number
        double ux, uy, uz; // Voxel spacing in real-world coordinates
        int nx, ny, nz, nc; // Image dimensions
        bool valid; // Data validity 

public:

        /* Data is initially invalid */
        Dicom() : filename(""), patient(-1), study(-1), series(-1), 
                instance(-1), valid(false) {};

        ~Dicom() {};

        /* Load a file */
        Dicom(std::string filename);

        /* Get the x-dimension */
        int getNx(void) const {
                return nx;
        }

        /* Get the y-dimension */
        int getNy(void) const {
                return ny;
        }

        /* Get the z-dimension */
        int getNz(void) const {
                return nz;
        }

        /* Get the number of channels */
        int getNc(void) const {
                return nc;
        } 

        /* Get the x-spacing */
        double getUx(void) const {
                return ux;
        }

        /* Get the y-spacing */
        double getUy(void) const {
                return uy;
        }

        /* Get the z-spacing */
        double getUz(void) const {
                return uz;
        }

        /* Check whether or not the data is valid */
        bool isValid(void) const {
                return valid;
        }

        /* Get the file name */
        std::string name(void) const {
                return filename;
        }

        /* Sort by instance number */
        bool operator < (const Dicom &dicom) const {
                return instance < dicom.instance;
        }

        /* Check if another DICOM file is from the same series */
        bool eqSeries(const Dicom &dicom) const {
                return patient == dicom.patient &&
                        study == dicom.study &&
                        series == dicom.series;
        }
};

/* Load the data from a DICOM file */
Dicom::Dicom(std::string path) : filename(path), valid(false) {

        // Load the image as a DcmFileFormat 
        DcmFileFormat fileFormat;
        OFCondition status = fileFormat.loadFile(path.c_str());
        if (!status.good()) {
               std::cerr << "Dicom.Dicom: failed to read DICOM file " <<
                        path << " (" << status.text() << ")" << 
                        std::endl; 
                return;
        }

        // Get the dataset
        DcmDataset *const data = fileFormat.getDataset();

        // Get the patient ID
        const char *patientIdStr;
        status = data->findAndGetString(DCM_PatientID, patientIdStr);
        if (status.bad()) {
                std::cerr << "Dicom.Dicom: failed to get patient ID " <<
                        "from file " << path << " (" << status.text() << ")" <<
                        std::endl;
                return;
        }
        patient = atoll(patientIdStr);

        // Get the study ID
        const char *studyIdStr;
        status = data->findAndGetString(DCM_StudyID, studyIdStr);
        if (status.bad()) {
                std::cerr << "Dicom.Dicom: failed to get study ID " <<
                        "from file " << path << " (" << status.text() << ")" <<
                        std::endl;
                return;
        }
        study = atoll(studyIdStr);

        // Get the series number
        const char *seriesStr;
        status = data->findAndGetString(DCM_SeriesNumber, seriesStr);
        if (status.bad()) {
                std::cerr << "Dicom.Dicom: failed to get series ID " <<
                        "from file " << path << " (" << status.text() << ")" <<
                        std::endl;
                return;
        }
        series = atoll(seriesStr); 

        // Get the instance number
        const char *instanceStr;
        status = data->findAndGetString(DCM_InstanceNumber, instanceStr);
        if (status.bad()) {
                std::cerr << "Dicom.Dicom: failed to get instance number " <<
                        "from file " << path << " (" << status.text() << ")" <<
                        std::endl;
                return;
        }
        instance = atoll(instanceStr);

        // Load the DicomImage object
        DicomImage image(path.c_str());
        if (image.getStatus() != EIS_Normal) {
               std::cerr << "Dicom.image: failed to open image " <<
                        filename << " (" << 
                        DicomImage::getString(image.getStatus()) << ")" << 
                        std::endl; 
                return;
        }

        // Check for color images
        if (!image.isMonochrome()) {
                std::cerr << "Dicom.Dicom: reading of color DICOM " <<
                        "images is not supported at this time" << std::endl;
                return;
        }
        nc = 1;

        // Read the dimensions
        nx = image.getWidth();
        ny = image.getHeight();
        nz = image.getFrameCount();
        if (nx < 1 || ny < 1 || nz < 1) {
                std::cerr << "Dicom.Dicom: invalid dimensions for file "
                        << filename << "(" << nx << ", " << ny << ", " << 
                        nz << ")" << std::endl;
                return;
        }

        // Read the pixel spacing
        Float64 pixelSpacing;
        status = data->findAndGetFloat64(DCM_PixelSpacing,
                pixelSpacing);
        if (status.bad()) {
                std::cerr << "Dicom.Dicom: failed to get pixel spacing " <<
                        "from file " << path << " (" << status.text() << ")" <<
                        std::endl;
                return;
        }
        ux = static_cast<double>(pixelSpacing);
        if (ux <= 0.0) {
                std::cerr << "Dicom.Dicom: file " << path << " has " <<
                        "invalid pixel spacing: " << ux << std::endl;
                return;
        }

        // Get the aspect ratio
        const double ratio = image.getHeightWidthRatio();
        uy = ux * ratio;
        if (uy <= 0.0) {
                std::cerr << "Dicom.Dicom: file " << path << " has invalid " <<
                        "pixel aspect ratio: " << ratio << std::endl;
                return;
        }

        // Read the slice thickness 
        Float64 sliceThickness;
        status = data->findAndGetFloat64(DCM_SliceThickness, sliceThickness);
        if (!status.good()) {
                std::cerr << "Dicom.Dicom: failed to get slice thickness " <<
                        "from file " << path << " (" << status.text() << ")" <<
                        std::endl;
                return;
        }

        // Convert to double 
        uz = sliceThickness;
        if (uz <= 0.0) {
                std::cerr << "Dicom.Dicom: file " << path << " has " <<
                        "invalid slice thickness: " << uz << std::endl;
                return;
        }
        
        // Set the window 
        image.setMinMaxWindow();

        valid = true;
}

/* Read a DICOM file into an Image struct. */
int read_dcm(const char *path, Image *const im) {

        int ret;

        CATCH_EXCEPTIONS(ret, "read_dcm", read_dcm_cpp, path, im);

        return ret;
}

/* Read all of the DICOM files from a directory into an Image struct. Slices 
 * must be ordered alphanumerically, starting with z = 0. */
int read_dcm_dir(const char *path, Image *const im) {

        int ret;

        CATCH_EXCEPTIONS(ret, "read_dcm_dir", read_dcm_dir_cpp, path, im);

        return ret;
}

/* Write an Image struct into a DICOM file. */
int write_dcm(const char *path, const Image *const im) {

        int ret;

        CATCH_EXCEPTIONS(ret, "write_dcm", write_dcm_cpp, path, im);

        return ret;
}

/* Write an Image struct into a directory of DICOM files. */
int write_dcm_dir(const char *path, const Image *const im) {

        int ret;

        CATCH_EXCEPTIONS(ret, "write_dcm_dir", write_dcm_dir_cpp, path, im);

        return ret;
}

/* Helper function to read a DICOM file using C++ */
static int read_dcm_cpp(const char *path, Image *const im) {

        // Read the image metadata
        Dicom dicom(path);
        if (!dicom.isValid())
                return SIFT3D_FAILURE;

        // Load the DicomImage object
        DicomImage image(path);
        if (image.getStatus() != EIS_Normal) {
               std::cerr << "read_dcm_cpp: failed to open image " <<
                        dicom.name() << " (" << 
                        DicomImage::getString(image.getStatus()) << ")" << 
                        std::endl; 
                return SIFT3D_FAILURE;
        }

        // Initialize the image fields
        im->nx = dicom.getNx();
        im->ny = dicom.getNy();
        im->nz = dicom.getNz();
        im->nc = dicom.getNc();
        im->ux = dicom.getUx();
        im->uy = dicom.getUy();
        im->uz = dicom.getUz();

        // Resize the output
        im_default_stride(im);
        if (im_resize(im))
                return SIFT3D_FAILURE;

        // Read each frame
        for (int i = 0; i < im->nz; i++) { 

                // Get a pointer to the data, rendered as a 32-bit int
                const uint32_t *const frameData = 
                        static_cast<const uint32_t *const >(
                                image.getOutputData(32, i));
                if (frameData == NULL) {
                        std::cerr << "read_dcm_dir_cpp: could not get data "
                                << "from image " << path << " frame " << i <<
                                " (" << 
                                DicomImage::getString(image.getStatus()) << 
                                ")" << std::endl; 
                        return SIFT3D_FAILURE;
                }

                // Copy the frame
                const int x_start = 0;
                const int y_start = 0;
                const int z_start = i;
                const int x_end = im->nx - 1;
                const int y_end = im->ny - 1;
                const int z_end = z_start;
                int x, y, z;
                SIFT3D_IM_LOOP_LIMITED_START(im, x, y, z, x_start, x_end, 
                        y_start, y_end, z_start, z_end)

                        SIFT3D_IM_GET_VOX(im, x, y, z, 0) =
                                static_cast<float>(frameData[x + y * im->nx]);

                SIFT3D_IM_LOOP_END
        }

        return SIFT3D_SUCCESS;
}

/* Helper funciton to read a directory of DICOM files using C++ */
static int read_dcm_dir_cpp(const char *path, Image *const im) {

        struct stat st;
        DIR *dir;
        struct dirent *ent;
        int i, nx, ny, nz, nc, num_files, off_z;

        // Verify that the directory exists
	if (stat(path, &st)) {
                std::cerr << "read_dcm_dir_cpp: cannot find file " << path <<
                        std::endl;
                return SIFT3D_FAILURE;
	} else if (!S_ISDIR(st.st_mode)) {
                std::cerr << "read_dcm_dir_cpp: file " << path <<
                        " is not a directory" << std::endl;
                return SIFT3D_FAILURE;
	}

        // Open the directory
        if ((dir = opendir(path)) == NULL) {
                std::cerr << "read_dcm_dir_cpp: unexpected error opening " <<
                        "directory" << std::endl;
                return SIFT3D_FAILURE;
        }

        // Get all of the .dcm files in the directory
        std::vector<Dicom> dicoms;
        while ((ent = readdir(dir)) != NULL) {

                const char *filename = ent->d_name;

                // Check if it is a DICOM file 
                if (im_get_format(filename) != DICOM)
                        continue;

                // Form the whole file path
                std::string fullfile = std::string(path) + std::string(&sep) +
                        filename;

                // Read the file
                Dicom dicom(fullfile);
                if (!dicom.isValid()) {
                        closedir(dir);
                        return SIFT3D_FAILURE;
                }

                // Add the file to the list
                dicoms.push_back(dicom);
        }

        // Release the directory
        closedir(dir);
        
        // Get the number of files
        num_files = dicoms.size();

        // Verify that dicom files were found
        if (num_files == 0) {
                std::cerr << "read_dcm_dir_cpp: no dicom files found in " <<
                        path << std::endl;
                return SIFT3D_FAILURE;
        }

        // Check that the files are from the same series
        const Dicom &first = dicoms[0];
        for (int i = 1; i < num_files; i++) {

                const Dicom &dicom = dicoms[i];

                if (!first.eqSeries(dicom)) {
                        std::cerr << "read_dcm_dir_cpp: file " << 
                                dicom.name() << 
                                "is from a different series than file " <<
                                first.name() << std::endl;
                        return SIFT3D_FAILURE;
                }
        }

        // Initialize the output dimensions
        nx = first.getNx();
        ny = first.getNy();
        nc = first.getNc();

        // Verify the dimensions of the other files, counting the total
        // series z-dimension
        nz = 0;
        for (i = 0; i < num_files; i++) {

                // Get a slice
                const Dicom &dicom = dicoms[i];        

                // Verify the dimensions
                if (dicom.getNx() != nx || dicom.getNy() != ny || 
                        dicom.getNc() != nc) {
                        std::cerr << "read_dcm_dir_cpp: slice " << 
                                dicom.name() <<
                                " (" << dicom.getNx() << "x, " << 
                                dicom.getNy() << "y, " << dicom.getNc() << 
                                "c) does not match the dimensions of slice " <<
                                first.name() << "(" << nx << "x, " << ny << 
                                "y, " << nc << "c). " << std::endl;
                        return SIFT3D_FAILURE;
                }

                // Count the z-dimension
                nz += dicom.getNz();
        }

        // Resize the output
        im->nx = nx;
        im->ny = ny;
        im->nz = nz;
        im->nc = nc;
        im->ux = first.getUx(); 
        im->uy = first.getUy();
        im->uz = first.getUz();
        im_default_stride(im);
        if (im_resize(im))
                return SIFT3D_FAILURE;

        // Sort the slices by instance number
        std::sort(dicoms.begin(), dicoms.end()); 

        // Allocate a temporary image for the slices
        Image slice;
        init_im(&slice);

        // Read the image data
        off_z = 0;
        for (i = 0; i < num_files; i++) {

                int x, y, z, c;

                const char *slicename = dicoms[i].name().c_str();

                // Read the slice 
                if (read_dcm(slicename, &slice)) {
                        im_free(&slice);
                        return SIFT3D_FAILURE;
                }

                // Copy the data to the volume
                SIFT3D_IM_LOOP_START_C(&slice, x, y, z, c)

                        SIFT3D_IM_GET_VOX(im, x, y, z + off_z, c) =
                                SIFT3D_IM_GET_VOX(&slice, x, y, z, c);

                SIFT3D_IM_LOOP_END_C

                off_z += slice.nz;
        }
        assert(off_z == nz);
        im_free(&slice);

        return SIFT3D_SUCCESS;
}

/* Helper function to write a DICOM file using C++ */
static int write_dcm_cpp(const char *path, const Image *const im) {

        // Ensure the image is monochromatic
        if (im->nc != 1) {
                std::cerr << "write_dcm_cpp: image has " << im->nc <<
                        " channels. Currently only single-channel images " <<
                        "are supported." << std::endl;
                return SIFT3D_FAILURE;
        }

        // Create a new fileformat object
        DcmFileFormat fileFormat;

        // Set the file type to derived
        DcmDataset *const dataset = fileFormat.getDataset();
        OFCondition status = dataset->putAndInsertString(DCM_ImageType, 
                                                         "DERIVED");
        if (status.bad()) {
                std::cerr << "write_dcm_cpp: Failed to set the image type" <<
                        std::endl;
                return SIFT3D_FAILURE;
        }

        // Convert the dimensions to strings
#define BUF_LEN 1024
        char buf[BUF_LEN];

        // Set the dimensions
        snprintf(buf, BUF_LEN, "%d", im->nx);
        OFCondition xstatus = dataset->putAndInsertString(DCM_Rows, buf); 
        snprintf(buf, BUF_LEN, "%d", im->ny);
        OFCondition ystatus = dataset->putAndInsertString(DCM_Columns, buf);
        snprintf(buf, BUF_LEN, "%d", im->nz);
        OFCondition zstatus = dataset->putAndInsertString(DCM_NumberOfFrames,
                buf);
        if (xstatus.bad() || ystatus.bad() || zstatus.bad()) {
                std::cerr << "write_dcm_cpp: Failed to set the dimensions " <<
                        std::endl;
                return SIFT3D_FAILURE;
        }

        // Set the pixel spacing
        snprintf(buf, BUF_LEN, "%f", im->ux);
        status = dataset->putAndInsertString(DCM_PixelSpacing, buf);
        if (status.bad()) {
                std::cerr << "write_dcm_cpp: Failed to set the pixel " <<
                        "spacing" << std::endl;
                return SIFT3D_FAILURE;
        }

        // Set the aspect ratio
        const double ratio = im->uy / im->ux;
        snprintf(buf, BUF_LEN, "%f", ratio);
        status = dataset->putAndInsertString(DCM_PixelAspectRatio, buf);
        if (status.bad()) {
                std::cerr << "write_dcm_cpp: Failed to set the pixel " <<
                        "aspect ratio" << std::endl;
                return SIFT3D_FAILURE;
        }

        // Set the slice thickness
        snprintf(buf, BUF_LEN, "%f", im->uz);
        status = dataset->putAndInsertString(DCM_SliceThickness, buf);
        if (status.bad()) {
                std::cerr << "write_dcm_cpp: Failed to set the slice " <<
                        "thickness" << std::endl;
                return SIFT3D_FAILURE;
        }

        // Count the number of pixels in the image
        unsigned long numPixels = 1;
        for (int i = 0; i < IM_NDIMS; i++) {
                numPixels *= im->dims[i];
        }

        // Render the data to an 8-bit unsigned integer array
        uint8_t *pixelData = new uint8_t[numPixels];
        int x, y, z;
        SIFT3D_IM_LOOP_START(im, x, y, z)
                pixelData[x + y * im->nx + z * im->nx * im->ny] =
                        static_cast<uint8_t>(
                        SIFT3D_IM_GET_VOX(im, x, y, z, 0));
        SIFT3D_IM_LOOP_END

        // Write the data
        status = dataset->putAndInsertUint8Array(DCM_PixelData, pixelData, 
                numPixels);
        delete pixelData;
        if (status.bad()) {
                std::cerr << "write_dcm_cpp: Failed to set the pixel data " <<
                        std::endl;
                return SIFT3D_FAILURE;
        }

        // Choose the encoding format
#if 0
        const E_TransferSyntax xfer = EXS_JPEGProcess14SV1TransferSyntax;
        DJ_RPLossless rp_lossless;
        status = dataset->chooseRepresentation(xfer, &rp_lossless);
#else
        const E_TransferSyntax xfer = EXS_LittleEndianExplicit;
        dataset->chooseRepresentation(xfer, NULL);
#endif
        if (!dataset->canWriteXfer(xfer)) {
                std::cerr << "write_dcm_cpp: Failed to choose the encoding " <<
                        "format " << std::endl;
                return SIFT3D_FAILURE;
        }

        // Save the file
        status = fileFormat.saveFile(path, xfer);
        if (status.bad()) {
                std::cerr << "write_dcm_cpp: Failed to write file " <<
                        path << " (" << status.text() << ")" << std::endl;
                return SIFT3D_FAILURE;
        }

        return SIFT3D_SUCCESS;
#undef BUF_LEN
}

/* Helper function to write an image to a directory of DICOM files using C++ */
static int write_dcm_dir_cpp(const char *path, const Image *const im) {

        //TODO copy each frame to a dummy image, write it with write_dcm

        return SIFT3D_SUCCESS;
}

