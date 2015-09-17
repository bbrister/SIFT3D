# SIFT3D Windows Installation Instructions

Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

## Installing from binaries

You can install SIFT3D in one of two ways, from binaries or from source. The easiest way is to install from binaries. Simply visit our [binary distributions](https://github.com/bbrister/SIFT3D-installers) repository, download the appropriate installer for your system, and run it.

## Installing from source

*This is a difficult process, recommended only for advanced users.*

This program has been successfully compiled and executed on the following Windows platforms:
* Windows 8 64-bit, using [TDM-GCC](http://tdm-gcc.tdragon.net/) 5.10 and CMake 3.3.1

Please follow the instructions below to compile and install SIFT3D from source.

1. Install [CMake](http://www.cmake.org).

2. Install MinGW via [TDM-GCC](http://tdm-gcc.tdragon.net/)
	a. Download the installer
	b. Run the installer and remember to select C, C++, and Fortran

3. Use [gnumex](http://gnumex.sourceforge.net/documentation.html#L131) to hack Matlab to compile MEX files with MinGW. *This needed only for the Matlab toolbox. Non-Matlab users can skip this step.*
	a. Download and extract gnumex
	b. Run Matlab as an administrator
	c. Run the "gnumex" program from within Matlab
	d. Make sure the required paths are detected
	e. Generate the files
	f. Run "mex -setup" from within Matlab. If all goes well, GCC will show up as an option. Follow the prompt to set up MEX with GCC.

4. Install [LAPACK for Windows](http://icl.cs.utk.edu/lapack-for-windows/lapack/index.html#libraries).
	1. Download the binaries for MinGW and your verison of Windows (lapack.dll, blas.dll)
	2. Move lapack.dll and blas.dll to the TDM-GCC/bin directory 

5. Install [zlib](http://zlib.net/).
	1. Download and extract the most recent version
	2. Use the CMake GUI to generate MinGW Makes
		1. Set the source folder to the location of zlib
		2. Generate -> MinGW Makefiles
	3. Compile and install with MinGW
		1. cd to the build directory
		2. "mingw32-make"
		3. "mingw32-make install" (may require administrator priveleges)

6. Install [nifticlib](http://sourceforge.net/projects/niftilib/files/nifticlib/).
	1. Download and extract the newest version
	2. Generate MinGW Makefiles with CMake
		1. Set the source folder to the location of nifticlib
		2. Generate -> MinGW Makefiles
	3. Compile and install with MinGW
		1. cd to build directory
		2. "mingw32-make"
		3. "mingw32-make install" (may require administrator priveleges)

7. Install SIFT3D.
	1. Download this repository
	2. Use the CMake GUI to generate MinGW Makes
		1. Navigate to SIFT3D
		2. Select "MinGW Makefiles"
		3. Configure. This will fail because CMake cannot find nifticlib.
		4. Set NIFTI_DIR to the location of nifticlib in your build
		5. Ensure that the Matlab libraries are .dll's and not .lib. You will need to manually edit the paths for Matlab_MEX_LIBRARY, Matlab_MX_LIBRARY, MWLAPACK_LIBRARY, and MWBLAS_LIBRARY, so that "libmex.lib" is changed to "libmex.dll", etc. This requires locating these files within your Matlab installation. Check the "bin" directories for .dll's.
		6. Generate -> MinGW Makefiles
	3. Compile and install with MinGW
		1. cd to the build directory
		2. "mingw32-make"
		3. "mingw32-make install" (may require administrator priveleges)

## Caveats

This program was originally developed for Unix-like platforms, so some features have been disabled in the Windows version. By default, the command line interface is not compiled on Windows systems. If you wish to compile it, you can set the CMake variable BUILD_CLI to ON. This is not officially supported, and the resulting executables may not function correctly.

The good news is that you can still access SIFT3D through the libraries and wrappers for other languages.
