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

In addition, if Matlab is detected on your system, the following Matlab functions are created:
- detectSift3D.m - Detect SIFT3D keypoints from a 3-dimensional array.
- extractSift3D.m - Extract SIFT3D descriptors from a 3-dimensional array.

See /examples for sample programs using the C and Matlab APIs.

## Dependencies

This code requires the following tools to compile:
- [CMake](http://www.cmake.org)
- A suitable C/C++ compiler, such as GCC or Clang/LLVM.

In addition, this code requires the following external libraries:
- [LAPACK](http://www.netlib.org/lapack/)
- [nifticlib](http://sourceforge.net/projects/niftilib/files/nifticlib/)

Please follow the instructions below to install the other dependencies for your specific system.

### Ubuntu Linux

As of version 14.04, the following command will install all dependencies and build tools:

	sudo apt-get install build-essential cmake liblapack-dev libnifti-dev

### Mac OSX

As of version 10.10.5, you can install the dependencies with [Homebrew](http://brew.sh/). First install Homebrew, if you haven't already, then run the following command:
 
        brew install cmake niftilib

## Installation instructions

On Unix-like systems, the following commands will generate Makefiles and use them to compile the binaries in a subdirectory called "build":

	mkdir build
	cd build
	cmake ..
	make

If for some reason CMake cannot find the dependencies, you can specify the paths manually with the Cmake GUI. 

Use the following command to install the files:

	sudo make install

In principle you can use CMake to compile this code on Windows, but some modifications may be required to resolve the external dependencies.

### Matlab toolbox 

If Matlab is detected in your system, a Matlab toolbox is compiled in the /build/wrappers/matlab subdirectory. To install the toolbox, add the following line to your startup.m file:

        run /path/to/SIFT3D/build/wrappers/matlab/setupSift3D

If you do not have a startup file, you can simply run this command prior to using any SIFT3D functions.

#### Known issues

CMake sometimes fails to find Matlab on your system, especially on Mac OSX. If that is the case, you should see "Matlab not found" after running the cmake command. You can fix this by manually specifying the location of Matlab in the variable Matlab_ROOT_DIR. For example,

        cmake .. -DMatlab_ROOT_DIR=/path/to/MATLAB

or on Mac OSX,

        cmake .. -DMatlab_ROOT_DIR=/path/to/Matlab.app

where /path/to/MATLAB is the location of your Matlab installation.

## Usage instructions

For instructions on using the CLI, use the "--help" option, e.g. 
        kpSift3D --help

### Linking to SIFT3D libraries with CMake

SIFT3D exports a CMake package to the install directories. Here is an example of linking to SIFT3D from a CMake list.

        find_package (SIFT3D) # Find SIFT3D
        add_executable (helloWorld helloWorld.c) # Declare a target
        target_link_libraries (helloWorld PUBLIC ${SIFT3D_LIBRARIES}) # Link to the SIFT3D libraries
        target_include_directories (helloWorld PUBLIC ${SIFT3D_INCLUDE_DIR}) # Find the SIFT3D headers

### Manually linking to SIFT3D libraries

Here is an example of compiling a C program by explicitly specifying the libraries and include directories:

```
gcc helloWorld.c -o helloWorld -I/usr/local/include/sift3d -L/usr/local/lib/sift3d -lreg -lsift3d -limutil -llapack -lblas -lz -lniftiio -lm
```

Linkage dependencies are as follows:
- libimutil - requires linking to zlib, nifticlib, LAPACK and BLAS
- libsift3d - requires linking to imutil
- libreg - requires linking to sift3d and imutil

## Contact

Please contact me at blaine@stanford.edu if you have any questions or concerns.

If you would like to cite this work, please refer to the following paper:

B. Rister, D. Reiter, H. Zhang, D. Volz, M. Horowitz, R. Gabr, and J. R. Cavallaro, "Scale- and Orientation-Invariant Keypoints in Higher-Dimensional Data," to appear in *Proceedings of the IEEE International Conference on Image Processing (ICIP)*. 2015.
