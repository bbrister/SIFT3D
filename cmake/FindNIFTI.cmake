# - FindNIFTI.cmake
#
# Author: Thomas Proeger
# Modified by Blaine Rister
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
FIND_PATH(NIFTI_INCLUDE_DIRS
NAMES nifti1.h
PATHS ${NIFTI_DIR}
HINTS "C:/Program Files (x86)/NIFTI" "C:/Program Files/NIFTI"
PATH_SUFFIXES include nifti include/nifti
DOC "The include directory containing nifti1.h"
)

# LIBRARY files
FIND_LIBRARY(NIFTI_NIFTICDF_LIBRARY
NAMES nifticdf
PATHS ${NIFTI_DIR}
HINTS "C:/Program Files (x86)/NIFTI" "C:/Program Files/NIFTI"
PATH_SUFFIXES lib
DOC "The library file libnifticdf.so"
)
FIND_LIBRARY(NIFTI_NIFTIIO_LIBRARY
NAMES niftiio
PATHS ${NIFTI_DIR}
HINTS "C:/Program Files (x86)/NIFTI" "C:/Program Files/NIFTI"
PATH_SUFFIXES lib
DOC "The library file libniftiiof.so"
)
FIND_LIBRARY(NIFTI_ZNZ_LIBRARY
NAMES znz
PATHS ${NIFTI_DIR}
HINTS "C:/Program Files (x86)/NIFTI" "C:/Program Files/NIFTI"
PATH_SUFFIXES lib
DOC "The library file libznz.so"
)

# Allow the user to specify the path to the NIFTI installation
get_filename_component (_NIFTI_DIR ${NIFTI_NIFTIIO_LIBRARY} DIRECTORY)
set (NIFTI_DIR ${_NIFTI_DIR} 
        CACHE PATH "The directory containing the NIFTI installation")

# The combined NIFTI libraries
set (NIFTI_LIBRARIES ${NIFTI_NIFTICDF_LIBRARY} ${NIFTI_NIFTIIO_LIBRARY} ${NIFTI_ZNZ_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set PNG_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NIFTI
"Cannot find package NIFTI. Did you install the package libnifti1-dev?"
NIFTI_INCLUDE_DIRS
NIFTI_NIFTICDF_LIBRARY
NIFTI_NIFTIIO_LIBRARY
NIFTI_ZNZ_LIBRARY
)
# these variables are only visible in 'advanced mode'
MARK_AS_ADVANCED(NIFTI_INCLUDE_DIRS
NIFTI_NIFTICDF_LIBRARY
NIFTI_NIFTIIO_LIBRARY
NIFTI_ZNZ_LIBRARY
)
