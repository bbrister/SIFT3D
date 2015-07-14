################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the wrappers for other programming languages.
################################################################################

# The name and location of the wrappers directory 
set (WRAPPERS_DIR "wrappers")
set (WRAPPERS_PATH "${CMAKE_BINARY_DIR}/${WRAPPERS_DIR}")

# Look for MATLAB
find_package (Matlab QUIET)
if (${Matlab_FOUND})
        message ("-- Matlab was found. Building the Matlab wrappers.")
endif ()

# Subdirectories for each language
if (${Matlab_FOUND})
        add_subdirectory (matlab)
endif ()
