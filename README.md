# SIFT3D
Analogue of the scale-invariant feature transform (SIFT) for three-dimensional images. Includes an image processing and linear algebra library with feature matching and RANSAC regression.

DEPENDENCIES

This code requires the following external libraries:
- LAPACK (http://www.netlib.org/lapack/)
- nifticlib (http://sourceforge.net/projects/niftilib/files/nifticlib/)

These will soon be replaced by ITK.
- ITK (http://www.itk.org/)

This code requires the following tools to compile:
- CMake (http://www.cmake.org)
- A suitable C/C++ compiler, such as GCC.

INSTALLATION INSTRUCTIONS

On Unix-like systems, the following commands will generate Makefiles and use them to compile the binaries in a subdirectory called "build."

mkdir build
cd build
cmake ..
make

In principle you can use CMake to compile this code on Windows, but some modifications may be required to resolve the external dependencies.

Please contact me at blaine@stanford.edu if you have any questions or concerns.
