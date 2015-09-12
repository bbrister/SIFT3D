################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Top-level build file for SIFT3D.
################################################################################

include (ExternalProject)

# Options
set (BUILD_EXAMPLES "ON" CACHE BOOL "If ON, builds the example programs.")
set (USE_Matlab "ON" CACHE BOOL "If ON, builds the Matlab toolbox")

# Minimum CMake version
cmake_minimum_required (VERSION 2.8)

# Project name
project (SIFT3D)

# Project version
set (SIFT3D_VERSION 1.0.0)

# Configurable paths        
set (INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib/sift3d" CACHE PATH 
        "Installation directory for libraries")
set (INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin/sift3d" CACHE PATH 
        "Installation directory for executables")
set (INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH 
        "Installation directory for header files")
if (WIN32 AND NOT CYGWIN)
        set (DEFAULT_INSTALL_CMAKE_DIR cmake)
else ()
        set (DEFAULT_INSTALL_CMAKE_DIR lib/cmake/sift3d)
endif ()
set (INSTALL_CMAKE_DIR ${DEFAULT_INSTALL_CMAKE_DIR} CACHE PATH
        "Installation directory for CMake files")

# Configure RPATH options
set (CMAKE_SKIP_BUILD_RPATH OFF)
set (CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
set (CMAKE_INSTALL_RPATH ${INSTALL_LIB_DIR})
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
LIST (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES ${INSTALL_LIB_DIR} 
        IS_SYSTEM_DIR)

# Internal paths 
list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
set (LICENSE_FILE ${CMAKE_CURRENT_LIST_DIR}/LICENSE)
set (README_FILE ${CMAKE_CURRENT_LIST_DIR}/README.md)
set (CMAKE_INSTALL_RPATH ${INSTALL_LIB_DIR})

# Output directories
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# The name and location of the wrappers build and install directories
set (WRAPPERS_DIR "wrappers")
set (BUILD_WRAPPERS_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${WRAPPERS_DIR}")
set (INSTALL_WRAPPERS_DIR "${INSTALL_LIB_DIR}/${WRAPPERS_DIR}")

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
if (WIN32)
        # m is not a separate library in Windows
        set (M_LIBRARY "")
else () 
        find_library (M_LIBRARY m)
endif ()
find_package (ZLIB REQUIRED)

# Optionally find Matlab
if (USE_Matlab)
	# Look for MATLAB in the default locations
	find_package (Matlab QUIET)

	# If not found, look for MATLAB in OS-specific locations
	if (NOT ${Matlab_FOUND})
		if (APPLE)
			set (Matlab_ROOT_DIR 
				"/Applications/Matlab.app")
			find_package (Matlab QUIET)
		endif ()
	endif ()
else ()
	# Disable Matlab
	set (Matlab_FOUND FALSE)
endif ()

# Set up Matlab paths
if (${Matlab_FOUND})
        message ("-- Found Matlab. Building the Matlab toolbox.")

        # Set the name and location of the Matlab toolbox
        set (TOOLBOX_NAME "matlab")
        set (BUILD_TOOLBOX_DIR "${BUILD_WRAPPERS_DIR}/${TOOLBOX_NAME}")
        set (INSTALL_TOOLBOX_DIR "${INSTALL_WRAPPERS_DIR}/${TOOLBOX_NAME}"
                CACHE PATH "Installation directory for the Matlab toolbox")

        # Get the path of the Matlab libraries
        get_filename_component (Matlab_LIBRARIES_PATH ${Matlab_MEX_LIBRARY} 
                DIRECTORY)

	# Check for the MX library, which CMake cannot find on Windows
	if (NOT Matlab_MX_LIBRARY_FOUND)

		# Find Matlab's MX library
		find_library (Matlab_MX_LIBRARY mx libmx 
                        PATHS ${Matlab_LIBRARIES_PATH}
                )
	endif ()

        # Find the Matlab LAPACK and BLAS libraries
        find_library (MWLAPACK_LIBRARY mwlapack libmwlapack 
                PATHS ${Matlab_LIBRARIES_PATH}
        )
        find_library (MWBLAS_LIBRARY mwblas libmwblas 
                PATHS ${Matlab_LIBRARIES_PATH}
        )
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

# CPack variables
set (CPACK_RESOURCE_FILE_LICENSE ${LICENSE_FILE})
set (CPACK_RESOURCE_FILE_README ${README_FILE})
set (CPACK_PACKAGE_VERSION ${SIFT3D_VERSION})
include (CPack)
