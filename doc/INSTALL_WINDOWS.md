# SIFT3D Windows Installation Instructions

Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.

## Installing from binaries

You can install SIFT3D in one of two ways, from binaries or from source. The easiest way is to install from binaries. Simply visit our [binary distributions](https://github.com/bbrister/SIFT3D-installers) repository, download the appropriate installer for your system, and run it.

## Installing from source

*This is a difficult process on Windows, recommended only for advanced users.*

This program has been successfully compiled and executed on the following Windows platforms:
* Windows 8 64-bit, using [TDM-GCC](http://tdm-gcc.tdragon.net/) 5.10 and CMake 3.3.1

Please follow the instructions below to compile and install SIFT3D from source.

1. Install [CMake](http://www.cmake.org).

2. Install MinGW via [TDM-GCC](http://tdm-gcc.tdragon.net/)
	1. Download the installer
	2. Run the installer and remember to select C, C++, and Fortran

3. Use [gnumex](http://gnumex.sourceforge.net/documentation.html#L131) to hack Matlab to compile MEX files with MinGW. *This needed only for the Matlab toolbox. Non-Matlab users can skip this step.*
	1. Download and extract gnumex
	2. Run Matlab as an administrator
	3. Run the "gnumex" program from within Matlab
	4. Make sure the required paths are detected
	5. Generate the files
	6. Run "mex -setup" from within Matlab. If all goes well, GCC will show up as an option. Follow the prompt to set up MEX with GCC.

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

7. Install [DCMTK](http://www.dcmtk.org).
	1. Download and extract the newest version. If there is a binary installer for your version of Windows, you can use that and skip the following sections. As of 10/06/2015, there is no binary for 64-bit Windows, so we compiled DCMTK from source. *Note: We could not compile DCMTK 3.6.0 on Windows with CMake and MinGW. Instead, we used a snapshot of version 3.6.1, retrieved on 09/24/2015.*
	2. Generate MinGW Makefiles with CMake
		1. Set the source folder to the location of DCMTK
		2. Generate -> MinGW Makefiles
	3. Compile and install with MinGW
		1. cd to build directory
		2. "mingw32-make"
		3. "mingw32-make install" (may require administrator priveleges)

8. Install SIFT3D.
	1. Download this repository
	2. Generate MinGW Makefiles with CMake
		1. Navigate to SIFT3D
		2. Select "MinGW Makefiles"
		3. Configure. This will fail because CMake cannot find DCMTK.
		4. Set DCMTK_DIR to the install location of DCMTK
		3. Configure again. This will fail because CMake cannot find nifticlib.
		4. Set NIFTI_DIR to the install location of nifticlib
		5. *Note: this step applies only to users compiling the optional Matlab toolbox.* Ensure that the Matlab libraries are .dll's and not .lib's. Manually edit the paths for Matlab_MEX_LIBRARY, Matlab_MX_LIBRARY, MWLAPACK_LIBRARY, and MWBLAS_LIBRARY, so that "libmex.lib" is changed to "libmex.dll", etc. This requires locating these files within your Matlab installation. Check the "bin" directories for .dll's.
		6. Generate -> MinGW Makefiles
	3. Compile and install with MinGW
		1. cd to the build directory
		2. "mingw32-make"
		3. "mingw32-make install" (may require administrator priveleges)

## Caveats

This program was originally developed for Unix-like platforms, so some features have been disabled in the Windows version. 

* By default, the command line interface is not compiled on Windows systems. If you wish to compile it, you can set the CMake variable BUILD_CLI to ON. This is not officially supported, and the resulting executables may not function correctly. The good news is that you can still access SIFT3D through the C libraries and wrappers for other languages.

* There are problems writing .nii.gz files in Windows. Note that all other image formats still work, including .nii without gzip compression. From the Matlab wrapper function imWrite3D, we work around this issue by first writing the .nii file, then compressing with a Matlab function call. However, you may experience issues when writing a .nii.gz file from the C library function im_write.
