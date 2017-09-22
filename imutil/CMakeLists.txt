################################################################################
# Copyright (c) 2015-2017 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the image processing utility library.
################################################################################

# Find BLAS
find_package (BLAS QUIET)
if (NOT BLAS_FOUND)
    message (FATAL_ERROR "BLAS not found. Please set the variable "
	    "BLAS_LIBRARIES to the location of the BLAS library on "
	    "your system.")
endif() 

# Find LAPACK
find_package (LAPACK QUIET)
if (NOT LAPACK_FOUND)
    message (FATAL_ERROR "LAPACK not found. Please set the variable "
	    "LAPACK_LIBRARIES to the location of the LAPACK library on your "
	    "system.")
endif ()

# Find DCMTK
find_package (DCMTK QUIET)
if (DCMTK_FOUND)

        # Add platform-specific DCMTK dependencies
        if (WIN32)
                # Add ws2_32 
                list (APPEND DCMTK_LIBRARIES "ws2_32")
        endif ()

        # Get the base DCMTK include dir. Note that DCMTK_DIR is set incorrectly on 
        # Linux, so we must add additional paths
        find_path (DCMTK_BASE_INCLUDE_PARENT_DIR "include/dcmtk"
                PATHS ${DCMTK_DIR} "${DCMTK_config_INCLUDE_DIR}/../../..")
        set (DCMTK_BASE_INCLUDE_DIR "${DCMTK_BASE_INCLUDE_PARENT_DIR}/include" 
                CACHE PATH "DCMTK include directory")

        if (_DCMTK_BASE_INCLUDE_PARENT_DIR STREQUAL
                "DCMTK_BASE_INCLUDE_PARENT_DIR-NOTFOUND")
                message (FATAL_ERROR "Failed to find the DCMTK include "
                        "directory. Please set the variable "
                        "DCMTK_BASE_INCLUDE_DIR to <DCMTK-INSTALL>/include, "
                        "or disable DICOM support by setting WITH_DICOM to "
                        "false.")
        endif ()

        # Add the base dir to the DCMTK include paths
        list(APPEND DCMTK_INCLUDE_DIRS ${DCMTK_BASE_INCLUDE_DIR})

        # Check if there is a configuration file for DCMTK
        find_file(DCMTK_CONFIG_FILE 
                NAMES "cfunix.h" "cfwin32.h"
                PATHS ${DCMTK_config_INCLUDE_DIR}
                NO_CMAKE_PATH
                NO_CMAKE_ENVIRONMENT_PATH
                NO_SYSTEM_ENVIRONMENT_PATH)
        if (DCMTK_CONFIG_FILE STREQUAL "DCMTK_CONFIG_FILE-NOTFOUND")
                set (DCMTK_HAVE_CONFIG_FILE false)
        else ()
                set (DCMTK_HAVE_CONFIG_FILE true)
        endif ()

	message (STATUS "Found DCMTK.")

elseif (WITH_DICOM)
	message (FATAL_ERROR "DCMTK not found. Please set the variable "
                "DCMTK_DIR to the location of DCMTK on your system, or disable"
                "DICOM support by setting WITH_DICOM to false.")
else ()
	message (STATUS "DCMTK not found. Compiling without DICOM support. "
                "To enable DICOM support, set the variable DCMTK_DIR to the "
                "location of DCMTK on your system.")
endif ()
set (WITH_DICOM ${DCMTK_FOUND} CACHE BOOL "Compile DICOM I/O support")
if (WITH_DICOM)
	message (STATUS "Building with DICOM support.")
endif ()

# Optionally find NIFTI 
find_package (NIFTI QUIET)
if (NIFTI_FOUND)
	message (STATUS "Found NIFTI.")
elseif (WITH_NIFTI)
	message (FATAL_ERROR "Failed to find nifticlib. Please set the "
                "variable NIFTI_DIR to the location of DCMTK on your system, "
                "or disable NIFTI support by setting WITH_NIFTI to false.")
else()
	message (STATUS "Failed to find nifticlib. Compiling without NIFTI "
                "support. To enable NIFTI support, set the variable NIFTI_DIR "
                "to the location of NIFTI on your system.")
endif()
set (WITH_NIFTI ${NIFTI_FOUND} CACHE BOOL "Compile NIFTI I/O support")
if (WITH_NIFTI)
	message (STATUS "Building with NIFTI support.")
endif ()

# Find iconv on Mac
if (APPLE)
	find_library (ICONV_LIBRARY NAMES iconv libiconv libiconv-2 c REQUIRED)
endif ()

# Find MinGW dependencies
if (MINGW)
        find_package (MINGW REQUIRED)
endif ()

# Macro to optionally add DCMTK to a target
macro (add_DICOM arg)
        if (WITH_DICOM)
                target_include_directories(${arg} PRIVATE ${DCMTK_INCLUDE_DIRS})
                target_link_libraries (${arg} PRIVATE ${DCMTK_LIBRARIES})
                target_compile_definitions(${arg} PRIVATE "SIFT3D_WITH_DICOM")
                if (DCMTK_HAVE_CONFIG_FILE)
                        target_compile_definitions (${arg} PRIVATE 
                                "HAVE_CONFIG_H")
                endif ()
        endif ()
endmacro ()

# Macro to optionally add NIFTI to a target
macro (add_NIFTI arg)
        if (WITH_NIFTI)
                target_compile_definitions(${arg} PRIVATE "SIFT3D_WITH_NIFTI")
                target_include_directories(${arg} PRIVATE ${NIFTI_INCLUDE_DIRS})
                target_link_libraries (${arg} PRIVATE ${NIFTI_LIBRARIES})
        endif ()
endmacro ()

# Format the compiler definitions
set (IMUTIL_DEFINITIONS "SIFT3D_VERSION_NUMBER=${SIFT3D_VERSION}")

# Check if the compiler uses strndup and strnlen
check_function_exists (strnlen HAVE_STRNLEN)
if (HAVE_STRNLEN)
        list (APPEND IMUTIL_DEFINITIONS "SIFT3D_HAVE_STRNLEN")
endif ()
check_function_exists (strndup HAVE_STRNDUP)
if (HAVE_STRNDUP)
        list (APPEND IMUTIL_DEFINITIONS "SIFT3D_HAVE_STRNDUP")
endif ()

# Compile imutil
add_library (imutil SHARED imutil.c nifti.c dicom.cpp)
target_include_directories (imutil PUBLIC 
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)
target_link_libraries (imutil PRIVATE ${LAPACK_LIBRARIES})
target_compile_definitions (imutil PRIVATE ${IMUTIL_DEFINITIONS})
add_DICOM(imutil)
add_NIFTI(imutil)
install (FILES imtypes.h immacros.h imutil.h kernels.cl
        DESTINATION ${INSTALL_INCLUDE_DIR})

# Link to system libraries
target_link_libraries(imutil PRIVATE ${ZLIB_LIBRARIES} ${M_LIBRARY})
target_include_directories(imutil PRIVATE ${ZLIB_INCLUDE_DIR})
if (APPLE)
	target_link_libraries(imutil PUBLIC ${ICONV_LIBRARY})
endif ()

# Configure the installation
install (TARGETS imutil 
        EXPORT SIFT3D-targets 
        RUNTIME DESTINATION ${INSTALL_BIN_DIR} 
	LIBRARY DESTINATION ${INSTALL_LIB_DIR} 
	ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
)

# OS-specific installation
if (WIN32)

	# Make a list of all external dependencies
	set (DEPS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${ZLIB_LIBRARIES} 
                ${OpenMP_C_LIBRARIES} ${MINGW_LIBRARIES})
        if (WITH_DICOM)
                list (APPEND DEPS ${DCMTK_LIBRARIES})
        endif ()
        if (WITH_NIFTI)
                list (APPEND DEPS ${NIFTI_LIBRARIES})
        endif ()

	function (get_runtime_deps ARG_DEPS ARG_RUNTIME_DEPS)

		# Process each dependency, adding a runtime dependency if 
		# necessary
		set (${ARG_RUNTIME_DEPS} "")
		foreach (DEP IN LISTS ${ARG_DEPS})

			# Get the file extension	
			get_filename_component (DEP_EXT ${DEP} EXT)		

			# Process shared libraries
			if (DEP_EXT STREQUAL ".dll")
				list (APPEND ${ARG_RUNTIME_DEPS} ${DEP})
			# Process MinGW import libraries
			elseif (DEP_EXT STREQUAL ".dll.a")

				# Extract the filename, parent and grandparent directories
				get_filename_component (DEP_NAME ${DEP} NAME)
				get_filename_component (DEP_DIR ${DEP} DIRECTORY)
				get_filename_component (DEP_DIR_DIR ${DEP_DIR} DIRECTORY)

				# Get the name of the .dll version
				string (REGEX REPLACE ".dll.a$" ".dll" DEP_DLL_NAME ${DEP_NAME})

				# Find the corresponding .dll
				string (REGEX REPLACE ".dll" "_DLL" DEP_DLL_VAR ${DEP_DLL_NAME})
				find_file (${DEP_DLL_VAR} ${DEP_DLL_NAME} 
					PATHS ${DEP_DIR} ${DEP_DIR_DIR}
					PATH_SUFFIXES "bin" "lib")
				if (${DEP_DLL_VAR} STREQUAL "${DEP_DLL_NAME}-NOTFOUND")
					message (FATAL_ERROR 
						"Failed to find runtime dependency ${DEP_DLL_NAME}")
				endif ()

				# The .dll, not the .dll.a, becomes a runtime dependency
				list (APPEND ${ARG_RUNTIME_DEPS} ${${DEP_DLL_VAR}})
			endif ()

		endforeach ()

		# Set the return value
		set (${ARG_RUNTIME_DEPS} ${${ARG_RUNTIME_DEPS}} PARENT_SCOPE)

	endfunction ()

	# Convert dependencies to runtime dependencies
	get_runtime_deps (DEPS RUNTIME_DEPS)

	# Copy the runtime dependencies to the Windows DLL
	file (COPY ${RUNTIME_DEPS} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

	# Add the runtime dependencies to the Windows installer
	install (FILES ${RUNTIME_DEPS} DESTINATION ${INSTALL_BIN_DIR})

endif ()

# If Matlab was found, compile a copy for use with Matlab libraries
if (BUILD_Matlab)

        add_library (meximutil SHARED imutil.c nifti.c dicom.cpp)
        target_compile_definitions (meximutil PUBLIC "SIFT3D_MEX")
        target_compile_definitions (meximutil PRIVATE  ${IMUTIL_DEFINITIONS})

        target_include_directories (meximutil PUBLIC 
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
        )
        target_include_directories(meximutil PUBLIC ${Matlab_INCLUDE_DIRS})
        target_include_directories (meximutil PRIVATE 
                ${ZLIB_INCLUDE_DIR})
        target_link_libraries (meximutil PUBLIC ${Matlab_LIBRARIES})
        target_link_libraries (meximutil PRIVATE ${Matlab_MWLAPACK_LIBRARY} 
                ${Matlab_MWBLAS_LIBRARY} ${ZLIB_LIBRARIES} ${M_LIBRARY})
        add_DICOM(meximutil)
        add_NIFTI(meximutil)

        set_target_properties (meximutil 
                PROPERTIES 
                ARCHIVE_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
                LIBRARY_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
                RUNTIME_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
        )

        install (TARGETS meximutil 
                RUNTIME DESTINATION ${INSTALL_TOOLBOX_DIR}
                LIBRARY DESTINATION ${INSTALL_TOOLBOX_DIR}
                ARCHIVE DESTINATION ${INSTALL_TOOLBOX_DIR}
        )

        if (WIN32)

		# The toolbox has the same dependencies except for BLAS and LAPACK
		set (TOOLBOX_DEPS ${DEPS})
		list (REMOVE_ITEM TOOLBOX_DEPS 
			${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
		)

		get_runtime_deps (TOOLBOX_DEPS TOOLBOX_RUNTIME_DEPS)
		file (COPY ${TOOLBOX_RUNTIME_DEPS} DESTINATION 
			${BUILD_TOOLBOX_DIR})
		install (FILES ${TOOLBOX_RUNTIME_DEPS} 
			DESTINATION ${INSTALL_TOOLBOX_DIR})
        endif ()
endif ()

# Add the code snippets
add_subdirectory (templates)
