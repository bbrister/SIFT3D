# SIFT3D Mac Installation Instructions

Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

# Installing the dependencies

This program requires the following external libraries:
- [LAPACK](http://www.netlib.org/lapack/)
- [nifticlib](http://sourceforge.net/projects/niftilib/files/nifticlib/)

As of version 10.10.5, you can install the dependencies with [Homebrew](http://brew.sh/). First install Homebrew, if you haven't already, then run the following command:
 
        brew install niftilib

# Installing SIFT3D

## Installing from binaries

You can install SIFT3D in one of two ways, from binaries or from source. The easiest way is to install from binaries. Simply visit our [binary distributions](x) *TODO* repository, download the appropriate installer for your system, and run it. 

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
	cmake ..
	make

If for some reason CMake cannot find the dependencies, you can specify the paths manually with the CMake GUI. 

Use the following command to install the files:

	sudo make install
