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

# Paths 
list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

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

set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${RELEASE_FLAGS}")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELEASE_FLAGS}")

# OS-specific build flags
if (APPLE)
        set (CMAKE_MACOSX_RPATH ON)
endif()

# Find the system libraries
find_library (M_LIBRARY m)
find_library (Z_LIBRARY z)

# Find LAPACK
find_package (LAPACK QUIET)

# Mac OSX workaround for LAPACK
if (APPLE)
        find_library (name LAPACK)
endif ()

if (NOT LAPACK_FOUND)
	message ("LAPACK was not found on your system. Attempting to install from source...")

	ExternalProject_Add (
		LAPACK	
		CMAKE_ARGS -DBUILD_SHARED_LIBS=ON
		GIT_REPOSITORY https://github.com/qsnake/lapack.git 
		INSTALL_COMMAND sudo make install
	)

	set (LAPACK_LIBRARIES "-llapack -lblas")
endif ()

# Try to find ITK
find_package (ITK 4 QUIET)

# Source subdirectories
add_subdirectory (imutil)
add_subdirectory (sift3d)
add_subdirectory (reg)
add_subdirectory (cli)
add_subdirectory (wrappers)
if (${BUILD_EXAMPLES})
        add_subdirectory (examples)
endif ()
