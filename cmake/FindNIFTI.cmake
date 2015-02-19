# - FindNIFTI.cmake
#
# Author: Thomas Proeger
#
# This cmake find module looks for the header files and libraries from the
# 'libnifti1-dev' package.
#
# The following variables will be exported:
#
# NIFTI_INCLUDE_DIR - the directory that contains nifti1.h
# NIFTI_NIFTICDF_LIBRARY - the libnifticdf library file
# NIFTI_NIFTIIO_LIBRARY - the libniftiio library file
# NIFTI_ZNZ_LIBRARY - the libznz library file
# NIFTI_FOUND - TRUE if and only if ALL other variables have
# correct values.
#
# INCLUDE directory
FIND_PATH(NIFTI_INCLUDE_DIR
NAMES nifti1.h
PATH_SUFFIXES nifti
DOC "The include directory containing nifti1.h"
)
# LIBRARY files
FIND_LIBRARY(NIFTI_NIFTICDF_LIBRARY
NAMES nifticdf
DOC "The library file libnifticdf.so"
)
FIND_LIBRARY(NIFTI_NIFTIIO_LIBRARY
NAMES niftiio
DOC "The library file libniftiiof.so"
)
FIND_LIBRARY(NIFTI_ZNZ_LIBRARY
NAMES znz
DOC "The library file libznz.so"
)
# handle the QUIETLY and REQUIRED arguments and set PNG_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NIFTI
"Cannot find package NIFTI. Did you install the package libnifti1-dev?"
NIFTI_INCLUDE_DIR
NIFTI_NIFTICDF_LIBRARY
NIFTI_NIFTIIO_LIBRARY
NIFTI_ZNZ_LIBRARY
)
# these variables are only visible in 'advanced mode'
MARK_AS_ADVANCED(NIFTI_INCLUDE_DIR
NIFTI_NIFTICDF_LIBRARY
NIFTI_NIFTIIO_LIBRARY
NIFTI_ZNZ_LIBRARY
)
