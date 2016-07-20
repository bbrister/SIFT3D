################################################################################
# SIFT3DPackage.cmake
################################################################################
# Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file to populate SIFT3D packages.
################################################################################

# Set the default package generator based on the OS 
if (WIN32)
        # NSIS for Windows
        set (_CPACK_GENERATOR "NSIS")
elseif (APPLE)
        # Bundle for Mac
        set (_CPACK_GENERATOR "BUNDLE")
elseif (UNIX)
        # Default to .tar.gz on Unix-like platforms
        set (_CPACK_GENERATOR "TGZ")

        # Override for specific Linux distributions
        if (CMAKE_SYSTEM_NAME MATCHES "Linux")
                # Try to read the Linux distribution
                if (EXISTS "/etc/issue")
                        file (READ "/etc/issue" LINUX_ISSUE)
                        # .deb for Ubuntu and Debian
                        if (${LINUX_ISSUE} MATCHES "Ubuntu" OR 
                                ${LINUX_ISSUE} MATCHES "Debian")
                                set (_CPACK_GENERATOR "DEB")
                        # .rpm for Fedora and OpenSuSE
                        elseif (${LINUX_ISSUE} MATCHES "Fedora" OR
                                ${LINUX_ISSUE} MATCHES "SUSE")
                                set (_CPACK_GENERATOR "RPM")
                        endif ()
                endif ()
        endif ()
else ()
        message (FATAL_ERROR "Unable to determine the default package generator for this operating system")
endif ()
set (CPACK_GENERATOR ${_CPACK_GENERATOR} CACHE STRING "The package generation 
        program")

# Global CPack variables
set (CPACK_PACKAGE_NAME "SIFT3D")
set (CPACK_PACKAGE_CONTACT "Blaine Rister blaine@stanford.edu")
set (CPACK_RESOURCE_FILE_LICENSE ${LICENSE_FILE})
set (CPACK_RESOURCE_FILE_README ${README_FILE})
set (CPACK_PACKAGE_VERSION ${SIFT3D_VERSION})
set (CPACK_PACKAGE_DESCRIPTION "Extracts keypoints and descriptors from 3D images. Also contians libraries for image processing and registration. Includes wrappers for Matlab.")

# Use the CMake install path, unless this is Windows, in which case this would
# break CPack
if (NOT WIN32)
        set (CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif ()

# Generator-specific CPack variables
if (CPACK_GENERATOR STREQUAL "NSIS")

        # Add the option to append SIFT3D to the path
        set (CPACK_NSIS_MODIFY_PATH ON)

elseif (CPACK_GENERATOR STREQUAL "BUNDLE")
        message(FATAL_ERROR "Bundle support not yet implemented")
        # TODO Do we need a bundle info file?
elseif (CPACK_GENERATOR STREQUAL "DEB")

        # Set the target platform in Debian terms
        if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
                set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64)
        elseif (CMAKE_SYSTEM_PROCESSOR STREQUAL "i686")
                set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE i386)
        else ()
                set (CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})
        endif () 
        message ("-- Configured Debian package for architecture ${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")

        # Set the standard Debian package dependencies
        set (CPACK_DEBIAN_PACKAGE_DEPENDS "libc6, zlib1g, libnifti2, libdcmtk5, liblapack3")

        # Set compiler-specific Debian package dependencies
        if (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR
               CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
                set (CPACK_DEBIAN_PACKAGE_DEPENDS 
                        "${CPACK_DEBIAN_PACKAGE_DEPENDS}, libgcc1")
        endif () 
elseif (CPACK_GENERATOR STREQUAL "RPM")
        message(FATAL_ERROR "RPM support not yet implemented")
        # TODO Add RPM package dependencies (see debian dependencies)
endif ()

# Finally include the CPack module
include (CPack)
