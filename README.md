# SIFT3D
Analogue of the scale-invariant feature transform (SIFT) for three-dimensional images. Includes an image processing and linear algebra library with feature matching and RANSAC regression.

CONTENTS

This code creates the following executables:
- kpSift3D - Extract keypoints and descriptors form a single image.
- regSift3D - Extract matches and a geometric transformation from two images. 

and the following libraries:
- libreg.so - Image registration from SIFT3D features
- libsift3d.so - Extract and match SIFT3D features
- libimutil.so - Utility library for image processing, regression and linear algebra.

Documentation can be found at:

See /examples for sample programs using the C library and CLI.

DEPENDENCIES

This code requires the following external libraries:
- LAPACK (http://www.netlib.org/lapack/)
- nifticlib (http://sourceforge.net/projects/niftilib/files/nifticlib/)

If the build system cannot find LAPACK, it will try to download and compile the source.

If you are using ITK, the build system will automatically find nifticlib there. Otherwise, you must install nifticlib yourself. (See the Ubuntu installation command below.)
- ITK (http://www.itk.org/)

This code requires the following tools to compile:
- CMake (http://www.cmake.org)
- A suitable C/C++ compiler, such as GCC or Clang/LLVM.

On Ubuntu, as of version 14.04, the following command will install all dependencies and build tools:

	sudo apt-get install build-essential cmake liblapack-dev libnifti-dev

INSTALLATION INSTRUCTIONS

On Unix-like systems, the following commands will generate Makefiles and use them to compile the binaries in a subdirectory called "build":

	mkdir build
	cd build
	cmake ..
	make

If for some reason CMake cannot find the dependencies, you can specify the paths manually with the cmake GUI. 

Use the following command to install the files:

	sudo make install

In principle you can use CMake to compile this code on Windows, but some modifications may be required to resolve the external dependencies.

Please contact me at blaine@stanford.edu if you have any questions or concerns.

USAGE INSTRUCTIONS

For instructions on using the CLI, use the "--help" option, e.g. 
        kpSift3D --help

Here is an example of compiling a C program with the libraries:

gcc helloWorld.c -o helloWorld -I/usr/local/include/sift3d -L/usr/local/lib/sift3d -lreg -lsift3d -limutil -llapack -lblas

Linkage dependencies are as follows:
	- libimutil - requires linking to lapack and BLAS
	- libsift3d - requires linking to libimutil
        - libreg - requires linking to libsift3d and libimutil
