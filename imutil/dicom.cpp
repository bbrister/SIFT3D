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
        Image im; // The image data
        std::string filename; // DICOM file name
        long long int patient; // Patient ID
        long long int study; // Study ID 
        long long int series; // Series number
        long long int instance; // Identifying instance number
        bool valid; // Data validity 

public:

        /* Data is initially invalid */
        Dicom() : filename(""), patient(-1), study(-1), series(-1), 
                instance(-1), valid(false) {
                init_im(&im);
        };

        Dicom(const Dicom &dicom) : filename(dicom.filename), 
                patient(dicom.patient), study(dicom.study), 
                series(dicom.series), instance(dicom.instance), 
                valid(dicom.valid) {

                init_im(&im);

                if (im_copy_data(dicom.getImage(), &im))
                        throw -1;
        }

        ~Dicom() {
                im_free(&im);
        };

        /* Load a file */
        Dicom(std::string filename);

        /* Get the data as a DicomImage object */
        const Image *getImage(void) const {
                return &im;
        }
        
        /* Get the x-dimension */
        int getNx(void) const {
                return im.nx;
        }

        /* Get the y-dimension */
        int getNy(void) const {
                return im.ny;
        }

        /* Get the z-dimension */
        int getNz(void) const {
                return im.nz;
        }

        /* Get the number of channels */
        int getNc(void) const {
                return im.nc;
        } 

        /* Get the x-spacing */
        double getUx(void) const {
                return im.ux;
        }

        /* Get the y-spacing */
        double getUy(void) const {
                return im.uy;
        }

        /* Get the z-spacing */
        double getUz(void) const {
                return im.uz;
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

        // Shallow copy of another Dicom's data--use with caution
        void shallowCopy(const Dicom &src) {
                im = src.im;
                filename = src.filename;
                patient = src.patient;
                study = src.study;
                series = src.series; 
                instance = src.instance;
                valid = src.valid;
        }

        // Release the memory for an image without freeing--use with caution
        void shallowRelease(void) {
                im.data = NULL;
        } 
};

/* Dont copy image data on a swap */
namespace std {
        template<>
        void swap(Dicom &d1, Dicom &d2) {

                Dicom temp; 

                // Shallow copies
                temp.shallowCopy(d1);
                d1.shallowCopy(d2);
                d2.shallowCopy(temp);

                // Prevent double-free
                temp.shallowRelease();
        }
}

/* Load the data from a DICOM file */
Dicom::Dicom(std::string path) : filename(path), valid(false) {

        init_im(&im);

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
        im.ux = static_cast<double>(pixelSpacing);
        if (im.ux <= 0.0) {
                std::cerr << "Dicom.Dicom: file " << path << " has " <<
                        "invalid pixel spacing: " << im.ux << std::endl;
                return;
        }
        im.uy = im.ux * image.getHeightWidthRatio();

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
        im.uz = sliceThickness;
        if (im.uz <= 0.0) {
                std::cerr << "Dicom.Dicom: file " << path << " has " <<
                        "invalid slice thickness: " << im.uz << std::endl;
                return;
        }
        
        // Set the window 
        image.setMinMaxWindow();

        // Resize the internal image data
        im.nx = image.getWidth();
        im.ny = image.getHeight();
        im.nz = image.getFrameCount();
        im.nc = 1;
        im_default_stride(&im);
        if (im_resize(&im))
                return;

        // Read each frame
        for (int i = 0; i < im.nz; i++) { 

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
                        return;
                }

                // Copy the frame
                const int x_start = 0;
                const int y_start = 0;
                const int z_start = i;
                const int x_end = im.nx - 1;
                const int y_end = im.ny - 1;
                const int z_end = z_start;
                int x, y, z;
                SIFT3D_IM_LOOP_LIMITED_START(&im, x, y, z, x_start, x_end, 
                        y_start, y_end, z_start, z_end)

                        SIFT3D_IM_GET_VOX(&im, x, y, z, 0) =
                                static_cast<float>(frameData[x + y * im.nx]);

                SIFT3D_IM_LOOP_END
        }

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

        // Read the image
        Dicom dicom(path);
        if (!dicom.isValid())
                return SIFT3D_FAILURE;

        // Copy the data
        return im_copy_data(dicom.getImage(), im);

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
        im_default_stride(im);
        if (im_resize(im))
                return SIFT3D_FAILURE;

        // Sort the slices by instance number
        std::sort(dicoms.begin(), dicoms.end()); 

        // Copy the data
        off_z = 0;
        for (i = 0; i < num_files; i++) {

                int x, y, z, c;

                const Image *const slice = dicoms[i].getImage();

                // Copy the data
                SIFT3D_IM_LOOP_START_C(slice, x, y, z, c)

                        SIFT3D_IM_GET_VOX(im, x, y, z + off_z, c) =
                                SIFT3D_IM_GET_VOX(slice, x, y, z, c);

                SIFT3D_IM_LOOP_END_C

                off_z += slice->nz;
        }
        assert(off_z == nz);

        return SIFT3D_SUCCESS;
}

/* Helper function to write a DICOM file using C++ */
static int write_dcm_cpp(const char *path, const Image *const im) {

        //TODO

        return SIFT3D_SUCCESS;
}

/* Helper function to write an image to a directory of DICOM files using C++ */
static int write_dcm_dir_cpp(const char *path, const Image *const im) {

        //TODO copy each frame to a dummy image, write it with write_dcm

        return SIFT3D_SUCCESS;
}

