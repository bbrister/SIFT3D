# SIFT3D Linux Installation Instructions

Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

# Installing the dependencies

This program requires the following external libraries:
- [LAPACK](http://www.netlib.org/lapack/)
- [nifticlib](http://sourceforge.net/projects/niftilib/files/nifticlib/)

On Ubuntu 14.04, the following command will install all dependencies and build tools:

	sudo apt-get install build-essential cmake liblapack-dev libnifti-dev

# Installing SIFT3D

## Installing from binaries

You can install SIFT3D in one of two ways, from binaries or from source. The easiest way is to install from binaries. Simply visit our [binary distributions](https://github.com/bbrister/SIFT3D-installers) repository, download the appropriate installer for your system, and run it. 

## Installing from source

This program has been successfully compiled and executed on the following Linux platforms:
- Ubuntu Linux 14.04, using GCC 4.8.4 and CMake 2.8.12.2.

This program requires the following tools to compile:
- [CMake](http://www.cmake.org)
- A suitable C/C++ compiler. such as GCC or Clang/LLVM.

On Ubuntu 14.04, the following command will install CMake and GCC:

        sudo apt-get install build-essential cmake

The following commands will generate Makefiles and use them to compile the binaries in a subdirectory called "build":

	mkdir build
	cd build
	cmake ..
	make

If for some reason CMake cannot find the dependencies, you can specify the paths manually with the CMake GUI. 

Use the following command to install the files:

	sudo make install
