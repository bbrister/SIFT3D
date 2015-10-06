# SIFT3D

Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

Analogue of the scale-invariant feature transform (SIFT) for three-dimensional images. Includes an image processing and linear algebra library with feature matching and RANSAC regression. Also includes IO functions supporting a variety of image formats.

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
        target_include_directories (helloWorld PUBLIC ${SIFT3D_INCLUDE_DIRS}) # Find the SIFT3D headers

### Manually linking to SIFT3D libraries

Here is an example of compiling a C program by explicitly specifying the libraries and include directories:

```
gcc helloWorld.c -o helloWorld -I/usr/local/include/sift3d -L/usr/local/lib/sift3d -lreg -lsift3d -limutil -llapack -lblas -lz -lnifticdf -lniftiio -lznz -lm
```

Linkage dependencies are as follows:
- libimutil - requires linking to zlib, nifticlib, LAPACK and BLAS
- libsift3d - requires linking to imutil
- libreg - requires linking to sift3d and imutil

## Contact

Please contact me at blaine@stanford.edu if you have any questions or concerns.

If you would like to cite this work, please refer to the following paper:

B. Rister, D. Reiter, H. Zhang, D. Volz, M. Horowitz, R. Gabr, and J. R. Cavallaro, "Scale- and Orientation-Invariant Keypoints in Higher-Dimensional Data," to appear in *Proceedings of the IEEE International Conference on Image Processing (ICIP)*. 2015.
