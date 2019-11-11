################################################################################
# Copyright (c) 2015-2017 Blaine Rister et al., see LICENSE for details.
################################################################################
# Top-level build file for SIFT3D.
################################################################################

include (ExternalProject)

# Minimum CMake version
cmake_minimum_required (VERSION 3.3)

# Minimum C++ version
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# Convert GNU archives to .lib files in Windows
if (WIN32)
	set (CMAKE_GNUtoMS ON
		CACHE BOOL 
		"Convert .dll.a files to .lib" 
	)
endif ()

# Project name
project (SIFT3D)

# Project version
set (SIFT3D_VERSION 1.4.6)

# Default build type
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message (STATUS  "Build type not specified. Defaulting to 'Release'.")
        set (CMAKE_BUILD_TYPE Release CACHE STRING "The type of build" FORCE)
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif ()

# Default settings for optional components
if (WIN32)
        set (_BUILD_CLI OFF)
else ()
        set (_BUILD_CLI ON)
endif()

# Optional component variables
set (BUILD_CLI ${_BUILD_CLI} CACHE BOOL 
        "If ON, builds the command line interface")
set (BUILD_EXAMPLES "ON" CACHE BOOL "If ON, builds the example programs")
set (BUILD_PACKAGE "OFF" CACHE BOOL "If ON, builds the package generator")

# Configurable paths        
set (INSTALL_LIB_DIR "lib/sift3d" CACHE PATH 
        "Installation directory for libraries")
set (INSTALL_BIN_DIR "bin" CACHE PATH 
        "Installation directory for executables")
set (INSTALL_INCLUDE_DIR "include/sift3d" CACHE PATH 
        "Installation directory for header files")
if (WIN32 AND NOT CYGWIN)
        set (DEFAULT_INSTALL_CMAKE_DIR "cmake")
else ()
        set (DEFAULT_INSTALL_CMAKE_DIR "lib/cmake/sift3d")
endif ()
set (INSTALL_CMAKE_DIR ${DEFAULT_INSTALL_CMAKE_DIR} CACHE PATH
        "Installation directory for CMake files")

# Internal paths 
list (APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
set (LICENSE_FILE ${CMAKE_CURRENT_LIST_DIR}/LICENSE)
set (README_FILE ${CMAKE_CURRENT_LIST_DIR}/README.md)

# Configure RPATH options
set (INSTALL_LIB_PATH "${CMAKE_INSTALL_PREFIX}/${INSTALL_LIB_DIR}")
set (CMAKE_SKIP_BUILD_RPATH OFF)
set (CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
set (CMAKE_INSTALL_RPATH ${INSTALL_LIB_PATH})
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

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
endif ()

# OS-specific build flags
if (APPLE)
        # Enable undefined shared library symbols
        set (CMAKE_SHARED_LINKER_FLAGS 
		"${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")

        # Use rpath
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

# Try to find and use OpenMP
find_package (OpenMP)
if (OPENMP_FOUND)
        # Nothing extra needs to be done here.
elseif (WITH_OpenMP)
        message (FATAL_ERROR "OpenMP not found. Please provide a compiler that "
                "supports OpenMP, or disable the OpenMP by setting "
                "the variable WITH_OpenMP to false.")
else ()
        message (STATUS "Failed to find OpenMP. Compiling without it.")
endif ()
set (WITH_OpenMP ${OPENMP_FOUND} CACHE BOOL "If ON, parallelizes with OpenMP "
        "in release mode")
if (WITH_OpenMP)
        set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
endif ()

# Look for MATLAB in the default locations
find_package (Matlab QUIET)

# If not found, look for MATLAB in OS-specific locations
if (NOT Matlab_FOUND)
        if (APPLE)
                set (Matlab_ROOT_DIR 
                        "/Applications/Matlab.app")
                find_package (Matlab QUIET)
        endif ()
endif ()

# Set up Matlab paths
if (Matlab_FOUND)
	message (STATUS "Found Matlab.")
elseif (BUILD_Matlab)
        message (FATAL_ERROR "Matlab not found. Please set the variable "
                "Matlab_ROOT_DIR to the location of your Matlab installation, "
                "or disable the Matlab toolbox by setting BUILD_Matlab to "
                "false.")
else ()
        message (STATUS "Matlab not found. The toolbox will not be built.")
endif ()
set (BUILD_Matlab ${Matlab_FOUND} CACHE BOOL "If ON, builds the Matlab toolbox")
if (BUILD_Matlab)
        message (STATUS "Building the Matlab toolbox.")

        # Set the name and location of the Matlab toolbox
        set (TOOLBOX_NAME "matlab")
        set (BUILD_TOOLBOX_DIR "${BUILD_WRAPPERS_DIR}/${TOOLBOX_NAME}")
        set (INSTALL_TOOLBOX_DIR "${INSTALL_WRAPPERS_DIR}/${TOOLBOX_NAME}"
                CACHE PATH "Installation directory for the Matlab toolbox")
        set (INSTALL_TOOLBOX_PATH "${CMAKE_INSTALL_PREFIX}/${INSTALL_TOOLBOX_DIR}")
        list (APPEND CMAKE_INSTALL_RPATH ${INSTALL_TOOLBOX_PATH})

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
        find_library (Matlab_MWLAPACK_LIBRARY mwlapack libmwlapack 
                PATHS ${Matlab_LIBRARIES_PATH}
        )
        find_library (Matlab_MWBLAS_LIBRARY mwblas libmwblas 
                PATHS ${Matlab_LIBRARIES_PATH}
        )
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

# Mandatory source subdirectories
add_subdirectory (imutil)
add_subdirectory (sift3d)
add_subdirectory (reg)
add_subdirectory (wrappers)

# Command line interface
if (BUILD_CLI)
        add_subdirectory (cli)
endif ()

# Examples
if (BUILD_EXAMPLES)
        add_subdirectory (examples)
endif ()

# Packager file
if (BUILD_PACKAGE)
        include (SIFT3DPackage)
endif ()
