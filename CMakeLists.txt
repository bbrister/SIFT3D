################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Top-level build file for SIFT3D.
################################################################################

include (ExternalProject)

# Options
set (BUILD_EXAMPLES "ON" CACHE BOOL "If ON, builds the example programs.")

# Minimum CMake version
cmake_minimum_required (VERSION 2.8)

# Project name
project (SIFT3D)

# Project version
set (SIFT3D_VERSION 1.0.0)

# Configurable paths        
set (INSTALL_LIB_DIR lib/sift3d CACHE PATH 
        "Installation directory for libraries")
set (INSTALL_BIN_DIR bin/sift3d CACHE PATH 
        "Installation directory for executables")
set (INSTALL_INCLUDE_DIR include CACHE PATH 
        "Installation directory for header files")
if (WIN32 AND NOT CYGWIN)
        set (DEFAULT_INSTALL_CMAKE_DIR cmake)
else ()
        set (DEFAULT_INSTALL_CMAKE_DIR lib/cmake/sift3d)
endif ()
set (INSTALL_CMAKE_DIR ${DEFAULT_INSTALL_CMAKE_DIR} CACHE PATH
        "Installation directory for CMake files")

# Internal paths 
list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# The name and location of the wrappers directory 
set (WRAPPERS_DIR "wrappers")
set (WRAPPERS_PATH "${CMAKE_BINARY_DIR}/${WRAPPERS_DIR}")

# The name and location of the Matlab toolbox
set (TOOLBOX_NAME "matlab")
set (TOOLBOX_PATH "${WRAPPERS_PATH}/${TOOLBOX_NAME}")

# Output directories
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Build flags
set (DEBUG_FLAGS "-g -DVERBOSE")
set (RELEASE_FLAGS "-O3 -DNDEBUG")

# GCC-specific flags
if (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR 
    CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (DEBUG_FLAGS "${DEBUG_FLAGS} -ggdb3")
    set (RELEASE_FLAGS "${RELEASE_FLAGS} -ffast-math")
endif ()

# OS-specific build flags
if (APPLE)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -undefined dynamic_lookup")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -undefined dynamic_lookup")
        set (CMAKE_MACOSX_RPATH ON)
endif ()

set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${RELEASE_FLAGS}")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELEASE_FLAGS}")

# Find the system libraries
find_library (M_LIBRARY m)
find_library (Z_LIBRARY z)

# Look for MATLAB in the default locations
find_package (Matlab QUIET)

# Look for MATLAB in OS-specific locations
if (APPLE AND NOT Matlab_ROOT_DIR)
        set (Matlab_ROOT_DIR "/Applications/Matlab.app")
        find_package (Matlab QUIET)
endif ()

if (${Matlab_FOUND})
        message ("-- Found Matlab. Building the Matlab toolbox.")
else ()
        message ("-- Matlab not found. The toolbox will not be built.")
endif ()

# Generate the package configuration files
set (CMAKE_CONFIG_FILE ${CMAKE_CURRENT_BINARY_DIR}/SIFT3DConfig.cmake)
set (CMAKE_VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/SIFT3DConfigVersion.cmake)
configure_file (SIFT3DConfig.cmake.in
        ${CMAKE_CONFIG_FILE} @ONLY
)
configure_file (SIFT3DConfigVersion.cmake.in
        ${CMAKE_VERSION_FILE} @ONLY)

# Install the package configuration files
install (FILES ${CMAKE_CONFIG_FILE} ${CMAKE_VERSION_FILE} 
        DESTINATION ${INSTALL_CMAKE_DIR})

# Install the targets file
install (EXPORT SIFT3D-targets DESTINATION ${INSTALL_CMAKE_DIR})

# Source subdirectories
add_subdirectory (imutil)
add_subdirectory (sift3d)
add_subdirectory (reg)
add_subdirectory (cli)
add_subdirectory (wrappers)
if (${BUILD_EXAMPLES})
        add_subdirectory (examples)
endif ()
