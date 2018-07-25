# SIFT3D Linux Installation Instructions

Copyright (c) 2015-2018 Blaine Rister et al., see LICENSE for details.

# Installing the dependencies

This program requires the following external libraries:
- [zlib](http://www.zlib.net/)
- [LAPACK](http://www.netlib.org/lapack/)

In addition, the following libraries add optional support for reading and writing DICOM and NIFTI files:
- [DCMTK](http://dicom.offis.de/dcmtk.php.en)
- [nifticlib](http://sourceforge.net/projects/niftilib/files/nifticlib/)

The optional I/O support will be compiled if DMCTK and nifticlib are found on your system. Otherwise, the programs will run as normal, throwing an error if you try to use an unsupported file format.

On Ubuntu 16.04, the following command will install all dependencies:

	sudo apt-get install zlib1g-dev liblapack-dev libdcmtk-dev libnifti-dev

# Installing SIFT3D

## Installing from binaries

You can install SIFT3D in one of two ways, from binaries or from source. The easiest way is to install from binaries. Simply visit our [releases](https://github.com/bbrister/SIFT3D/releases) page, download the appropriate installer for your system, and run it. 

## Installing from source

This program has been successfully compiled and executed on the following Linux platforms:
- Ubuntu Linux 16.04, using GCC 5.4.0 and CMake 3.5.1.

This program requires the following tools to compile:
- [CMake](http://www.cmake.org)
- A suitable C/C++ compiler, such as GCC or Clang/LLVM.

On Ubuntu 16.04, the following command will install CMake and GCC:

        sudo apt-get install build-essential cmake

The following commands will generate Makefiles and use them to compile the binaries in a subdirectory called "build":

        cd /path/to/SIFT3D # Go to the directory where you downloaded SIFT3D 
	mkdir build # Make a directory to store the binaries
	cd build
	cmake .. # Create the Makefiles
	make # Compile the program

If for some reason CMake cannot find the dependencies, you can specify the paths manually with the CMake GUI. 

Use the following command to install the files:

	sudo make install

### Packaging

To create your own binary package, invoke CMake with the BUILD_PACKAGE variable set to ON. Then build the 'package' target using the Makefiles. For example,

        cmake .. -DBUILD_PACKAGE=ON
        make package

### Troubleshooting

If you are using an older version of Ubuntu, then you will need to update [DCMTK](http://dicom.offis.de/dcmtk.php.en) to the latest version before compiling SIFT3D.

