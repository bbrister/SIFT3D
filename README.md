# SIFT3D

Copyright (c) 2015-2019 Blaine Rister et al., see LICENSE for details.

SIFT3D is an analogue of the scale-invariant feature transform (SIFT) for three-dimensional images. It leverages volumetric data and real-world units to detect keypoints and extract a robust description of their content. It can also perform 3D image registration by matching SIFT3D features and fitting geometric transformations with the RANSAC algorithm. All of this is implemented in a cross-platform C library, with wrappers for Matlab.

SIFT3D includes imutil, a utility library for image processing and linear algebra. This library performs file IO in a variety of medical imaging formats, including DICOM and NIFTI.

## Contents

This code creates the following executables:
- kpSift3D - Extract keypoints and descriptors from a single image.
- regSift3D - Extract matches and a geometric transformation from two images. 

and the following libraries:
- libreg.so - Image registration from SIFT3D features
- libsift3d.so - Extract and match SIFT3D features
- libimutil.so - Utility library for image processing, regression and linear algebra. Includes IO functions for DICOM and NIFTI file formats.

It also contains a Matlab toolbox for calling the library functions from Matlab scripts. See the README in /wrappers/matlab for more information.

## Installation instructions

See doc/INSTALL_\<PLATFORM\>.md for instructions on installing SIFT3D for your specific platform.

## Usage instructions

For instructions on using the CLI, use the "--help" option, e.g. 
        kpSift3D --help

See /examples for sample programs using the C and Matlab APIs.

The following sections describe how to link a program to the SIFT3D libraries.

### Linking to SIFT3D libraries with CMake

SIFT3D exports a CMake package to the install directories. Here is an example of compiling a C program with SIFT3D from a CMake list.

        find_package (SIFT3D) # Find SIFT3D
        add_executable (helloWorld helloWorld.c) # Declare a target
        target_link_libraries (helloWorld PUBLIC ${SIFT3D_LIBRARIES}) # Link to the SIFT3D libraries
        if (WIN32) # Find the SIFT3D headers
            target_include_directories (helloWorld PUBLIC "${SIFT3D_DIR}/../${SIFT3D_INCLUDE_DIRS}") 
        else()
            target_include_directories (helloWorld PUBLIC ${SIFT3D_INCLUDE_DIRS}) 
        endif()

### Linking to SIFT3D libraries without CMake

The header files and libraries are installed to "sift3d" subdirectories in your installation tree. On most systems, you will need to add these subdirectories to your include and linker search paths. You will also need to link to the dependencies listed below.

- libimutil - requires linking to LAPACK, BLAS, and zlib. Linking to DCMTK and nifticlib are optional.
- libsift3d - requires linking to imutil
- libreg - requires linking to sift3d and imutil

Information about the dependencies can be found in the installation instructions.

*Note: On Windows systems, some of the dependencies are statically linked to the SIFT3D libraries. In this case, it suffices to link to the DLLs in the "bin" subdirectory of your installation.*

## Contact

Please contact me at blaine@stanford.edu if you have any questions or concerns.

If you would like to cite this work, please refer to the following paper:

B. Rister, M. A. Horowitz and D. L. Rubin, "Volumetric Image Registration From Invariant Keypoints," in *IEEE Transactions on Image Processing*, vol. 26, no. 10, pp. 4900-4910, Oct. 2017.
doi: 10.1109/TIP.2017.2722689

The paper and automatic citations are available [here](http://ieeexplore.ieee.org/document/7967757/citations).
