#######################################################################
# Top-level build file for 3D SIFT
#######################################################################

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
set (DEBUG_FLAGS "-g -ggdb3")
set (RELEASE_FLAGS "-O3 -ffast-math")
if (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR 
    CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
    set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_DEBUG} ${RELEASE_FLAGS}")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_DEBUG} ${RELEASE_FLAGS}")
endif ()

# Find LAPACK
find_package (LAPACK REQUIRED)

# Try to find ITK
find_package (ITK 4 QUIET)

# Source subdirectories
add_subdirectory (imutil)
add_subdirectory (sift3d)
add_subdirectory (reg)
add_subdirectory (itkwrap)
