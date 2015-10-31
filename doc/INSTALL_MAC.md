# SIFT3D Mac Installation Instructions

Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

# Installing the dependencies

This program requires the following external libraries:
- [LAPACK](http://www.netlib.org/lapack/)
- [DCMTK](http://dicom.offis.de/dcmtk.php.en)
- [nifticlib](http://sourceforge.net/projects/niftilib/files/nifticlib/)

As of version 10.10.5, you can install the dependencies with [Homebrew](http://brew.sh/). First install Homebrew, if you haven't already, then run the following command:
 
        brew install dcmtk niftilib

# Installing SIFT3D

## Installing from binaries

We do not currently have binaries for SIFT3D on Mac. You must install from source.

## Installing from source

This program has been successfully compiled and executed on the following Mac platforms:
- Mac OSX 10.10.5, using Clang 6.1.0, LLVM 3.6.0 and CMake 3.3.1

This program requires the following tools to compile:
- [CMake](http://www.cmake.org)
- A suitable C/C++ compiler. such as GCC or Clang/LLVM.

Clang/LLVM comes pre-installed. Using [Homebrew](http://brew.sh), the following command will install CMake.

        brew install cmake

The following commands will generate Makefiles and use them to compile the binaries in a subdirectory called "build":

	mkdir build
	cd build
	cmake .. -DDCMTK_DIR=/usr/local
	make

*Note: If you installed DCMTK to some place other than /usr/local, set DCMTK_DIR to that path in the above CMake command.*

If CMake cannot find the dependencies, you can specify the paths manually with the CMake GUI.

Use the following command to install the files:

	sudo make install

