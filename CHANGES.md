# SIFT3D version history

## 1.0.0

* Initial release

## 1.1.0

* Added DICOM IO
* Wrapped DICOM and NIFTI IO in C functions im_read and im_write
* Added matlab wrappers keypoint3D, imRead3D and imWrite3D
* Changed C function SIFT3D_Extract_Descriptors to SIFT3D_Extract_Descriptors and SIFT3D_Extract_Raw_Descriptors
* Fixed various bugs, improving keypoint and descriptor accuracy

## 1.1.1 October 30, 2015

* Performance optimizations
* Fixes to DICOM slice ordering (im_read, imRead3D)
* Write more DICOM metadata (im_write, imWrite3D)
* Corrected Mac build instructions

## 1.2.0 January 5, 2016

* Fixed Mac linking issues
* Take real-world units into account in SIFT3D keypoints, descriptors
* Added Matlab wrapper for image registration
* Optionally resample input images prior to registration (regSift3D, register_SIFT3D_resample)
* Write more NIFTI metadata (im_write, imWrite3D)
* Changed C interface for descriptor extraction to disallow providing a custom Gaussian scale-space pyramid (SIFT3D_Extract_Descriptors)
* Allow arbitrary output sizes in im_inv_transform
* Minor bug fixes

## 1.3.0 March 28, 2016

* Add parameters for keypoint detection to Matlab interface
* Add registration parameters to Matlab and CLI
* Default in Matlab to faster version of image registration that avoids resampling
* Improved error handling
* Removed unused, faulty options
* Print internal error messages to the Matlab command prompt

## 1.3.1 March 31, 2016

* Removed option to set the number of octaves, which was causing bugs 

## 1.4.0 May 11, 2016

* Fixed bugs to improve keypoint and registration accuracy, especially concerning rotations
* Fixed bug that caused crashing on some large images
* Dramatically reduced memory consumption of keypoint matching
* Refactored SIFT3D_nn_match_fb into SIFT3D_nn_match, as there is now no reason to prefer forward to forward-backward matching
* Renamed headers macros.h and types.h to immacros.h and imtypes.h, respectively

## 1.4.1 May 25, 2016

* Removed keypoint refinement step, which did not improve the accuracy.

## 1.4.2 June 15, 2016

* Added multithreading with OpenMP
* Improved keypoint matching speed

## 1.4.3 July 24, 2016

* Fixed CMake settings to allow linking to CMake targets without dependencies
* Updated for newer versions of CMake, MinGW
* Update package dependencies for Ubuntu 16.04
* Removed POSIX dependencies in header files to allow linking with Visual Studio
* Added support for reading JPEG-compressed DICOM files
* Ship both MS (.lib) and MinGW (.dll.a) import libraries on Windows
* Ship with MinGW runtime libraries on Windows
* Ship with OpenBLAS on Windows

## 1.4.4 September 13, 2017

* Add support for reading Dicom Segmentation Objects (DSOs)
* Add the option to compile without DCMTK and nifticlib
* Reading images no longer scales them (im_read, imRead3D)
* Read DICOM CT scans in Hounsfield units (im_read, imRead3D)
* Fix header includes for newer builds of MinGW (TDM-GCC)

## 1.4.5

* Add wrapper function to match pre-computed descriptors from Matlab (matchSift3D)
