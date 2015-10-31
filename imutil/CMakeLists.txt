################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the image processing utility library.
################################################################################

# Find LAPACK
find_package (LAPACK REQUIRED)

# Find DCMTK
find_package (DCMTK 
	REQUIRED
	PATHS "/usr" "/usr/local" "/opt" "/opt/local" 
	PATH_SUFFIXES "include/dcmtk/config"
)

# Find NIFTI 
find_package (NIFTI REQUIRED)

# Try to find the path for DCMTK libraries
if (DCMTK_dcmdata_LIBRARY)
	get_filename_component (DCMTK_LIBRARIES_PATH 
		${DCMTK_dcmdata_LIBRARY} DIRECTORY)
endif ()

# Check if CMake is missing the new DCMTK library liboflog
if (NOT DCMTK_oflog_LIBRARY AND NOT TARGET oflog)

	# Search for liboflog
	find_library (DCMTK_oflog_LIBRARY oflog liboflog
		PATHS ${DCMTK_LIBRARIES_PATH})
	if (DCMTK_oflog_LIBRARY STREQUAL "DCMTK_oflog_LIBRARY-NOTFOUND")
		message (FATAL_ERROR "Could not find DCMTK library liboflog")
	endif ()

	# Find libofstd in the libraries list, which should be the last DCMTK 
	# library
	if (DEFINED DCMTK_LIBRARIES AND NOT TARGET ofstd)
		list (FIND DCMTK_LIBRARIES ${DCMTK_ofstd_LIBRARY} OFSTD_IDX)
		if (${OFSTD_IDX} LESS 0)
			message(FATAL_ERROR "Failed to find libofstd in 
				DCMTK_LIBRARIES")
		endif ()
	else ()
		set (OFSTD_IDX 0)
	endif ()

	# Insert liboflog just before libofstd
	list (INSERT DCMTK_LIBRARIES ${OFSTD_IDX} ${DCMTK_oflog_LIBRARY})
endif ()

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
	message (FATAL_ERROR "Failed to find the DCMTK include directory. Please set the variable DCMTK_BASE_INCLUDE_DIR to <DCMTK-INSTALL>/include")
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

# Compile imutil
add_library (imutil SHARED imutil.c dicom.cpp)
target_include_directories (imutil PUBLIC 
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)
target_link_libraries (imutil PUBLIC ${LAPACK_LIBRARIES} ${DCMTK_LIBRARIES} 
        ${NIFTI_LIBRARIES})
target_include_directories(imutil PRIVATE ${DCMTK_INCLUDE_DIRS} 
        ${NIFTI_INCLUDE_DIRS})
target_compile_definitions (imutil PRIVATE
        "SIFT3D_VERSION_NUMBER=${SIFT3D_VERSION}")
if (DCMTK_HAVE_CONFIG_FILE)
        target_compile_definitions (imutil PRIVATE "HAVE_CONFIG_H")
endif ()
install (FILES types.h macros.h imutil.h kernels.cl
        DESTINATION ${INSTALL_INCLUDE_DIR})

# Link to system libraries
target_link_libraries(imutil PUBLIC ${ZLIB_LIBRARIES} ${M_LIBRARY})
target_include_directories(imutil PRIVATE ${ZLIB_INCLUDE_DIRS})

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
	set (DEPS ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${DCMTK_LIBRARIES} 
                ${NIFTI_LIBRARIES} ${ZLIB_LIBRARIES})

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
if (${Matlab_FOUND})

        add_library (meximutil SHARED imutil.c dicom.cpp)
        target_compile_definitions (meximutil PUBLIC "SIFT3D_MEX")
        target_compile_definitions (meximutil PRIVATE 
                "SIFT3D_VERSION_NUMBER=${SIFT3D_VERSION}")
        if (DCMTK_HAVE_CONFIG_FILE)
                target_compile_definitions (meximutil PRIVATE "HAVE_CONFIG_H")
        endif ()

        target_include_directories (meximutil PUBLIC 
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
        )
        target_include_directories (meximutil PRIVATE 
                ${DCMTK_INCLUDE_DIRS} ${NIFTI_INCLUDE_DIRS} 
                ${Matlab_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
        target_link_libraries (meximutil PUBLIC ${Matlab_MWLAPACK_LIBRARY} 
                ${Matlab_MWBLAS_LIBRARY} ${Matlab_LIBRARIES} ${DCMTK_LIBRARIES}
                ${NIFTI_LIBRARIES} ${ZLIB_LIBRARIES} ${M_LIBRARY})

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

		set (TOOLBOX_DEPS ${DCMTK_LIBRARIES} ${NIFTI_LIBRARIES} 
                        ${ZLIB_LIBRARIES})
		get_runtime_deps (TOOLBOX_DEPS TOOLBOX_RUNTIME_DEPS)
		file (COPY ${TOOLBOX_RUNTIME_DEPS} DESTINATION 
			${BUILD_TOOLBOX_DIR})
                install (FILES ${TOOLBOX_RUNTIME_DEPS} 
			DESTINATION ${INSTALL_TOOLBOX_DIR})
        endif ()
endif ()

# Add the code snippets
add_subdirectory (templates)
