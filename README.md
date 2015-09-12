# SIFT3D

Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

Analogue of the scale-invariant feature transform (SIFT) for three-dimensional images. Includes an image processing and linear algebra library with feature matching and RANSAC regression.

## Contents

This code creates the following executables:
- kpSift3D - Extract keypoints and descriptors from a single image.
- regSift3D - Extract matches and a geometric transformation from two images. 

and the following libraries:
- libreg.so - Image registration from SIFT3D features
- libsift3d.so - Extract and match SIFT3D features
- libimutil.so - Utility library for image processing, regression and linear algebra.

In addition, if Matlab is detected on your system, a Matlab toolbox is automatically generated. See the README in /wrappers/matlab for more information.

## Installation instructions

### Supported platforms

This program has been successfully compiled and executed on the following platforms:
- Ubuntu Linux 14.04, using GCC 4.8.4 and CMake 2.8.12.2.
- Mac OSX 10.10.5, using Clang 6.1.0, LLVM 3.6.0 and CMake 3.3.1

Windows support is currently in the experimental stage. We plan to release binaries and installation instructions in the near future.

### Installing the dependencies

This code requires the following tools to compile:
- [CMake](http://www.cmake.org)
- A suitable C/C++ compiler, such as GCC or Clang/LLVM.

In addition, this code requires the following external libraries:
- [LAPACK](http://www.netlib.org/lapack/)
- [nifticlib](http://sourceforge.net/projects/niftilib/files/nifticlib/)

Please follow the instructions below to install the dependencies for your specific system.

#### Ubuntu Linux

As of version 14.04, the following command will install all dependencies and build tools:

	sudo apt-get install build-essential cmake liblapack-dev libnifti-dev

#### Mac OSX

As of version 10.10.5, you can install the dependencies with [Homebrew](http://brew.sh/). First install Homebrew, if you haven't already, then run the following command:
 
        brew install cmake niftilib

### Instaling SIFT3D 

On Unix-like systems, the following commands will generate Makefiles and use them to compile the binaries in a subdirectory called "build":

	mkdir build
	cd build
	cmake ..
	make

If for some reason CMake cannot find the dependencies, you can specify the paths manually with the CMake GUI. 

Use the following command to install the files:

	sudo make install

### Advanced installation

The following tutorials suggest advanced ways to install the software. They may be useful to some users, but are not required.

#### Out-of-source build

You do not need to compile SIFT3D in a subdirectory of the source tree. Instead, you can compile it in an arbitrary directory with the following commands:

        cmake /path/to/SIFT3D
        make

where /path/to/SIFT3D is the path to the SIFT3D source code.

#### Installation without root access

On some systems, you may not have access to install programs to protected locations. In this case, you can install SIFT3D to a different location by setting the CMAKE_INSTALL_PREFIX variable. For example,

        cmake /path/to/SIFT3D -DCMAKE_INSTALL_PREFIX=/path/to/install
        make
        make install

where /path/to/install is the path to the desired install directory. Note that other programs will not automatically find SIFT3D in this location.

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
