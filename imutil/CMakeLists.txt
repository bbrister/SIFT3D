################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the image processing utility library.
################################################################################

# Find LAPACK
find_package (LAPACK REQUIRED)

# Optionally disable use of ITK
set (USE_ITK OFF CACHE PATH 
        "Look for dependencies in the ITK installation, if one is found")

# Compile imutil
add_library (imutil SHARED imutil.c)
target_include_directories (imutil PUBLIC 
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)
target_link_libraries (imutil PUBLIC ${LAPACK_LIBRARIES})
install (FILES types.h macros.h imutil.h kernels.cl 
        DESTINATION ${INSTALL_INCLUDE_DIR})

# Try to find ITK, unless disabled
if (${USE_ITK})
        find_package (ITK 4 QUIET)
else ()
        set (ITK_FOUND OFF)
endif ()

if (${ITK_FOUND})
	message(STATUS "Linking to ITK")

	target_link_libraries(imutil PUBLIC ${ITK_LIBRARIES})
	target_include_directories(imutil PUBLIC ${ITK_INCLUDE_DIRS})

# Otherwise link directly to libniftiio
else ()
        if (${USE_ITK})
                message(STATUS "ITK not found. Linking to niftiio directly")
        endif ()

	find_package (NIFTI REQUIRED)

	target_link_libraries (imutil PUBLIC ${NIFTI_NIFTIIO_LIBRARY})
	target_include_directories (imutil PUBLIC ${NIFTI_INCLUDE_DIR})

endif ()

# Link to system libraries
target_link_libraries(imutil PUBLIC ${Z_LIBRARY} ${M_LIBRARY})

# If Matlab was found, compile a copy for use with Matlab libraries
if (${Matlab_FOUND})

        add_library (meximutil SHARED imutil.c)
        target_compile_definitions (meximutil PUBLIC "-DSIFT3D_MEX")
        set_target_properties (meximutil 
                PROPERTIES 
                ARCHIVE_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
                LIBRARY_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
                RUNTIME_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        )
        target_include_directories (meximutil PUBLIC 
                ${Matlab_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
        target_link_libraries (meximutil PUBLIC ${Z_LIBRARY} ${M_LIBRARY})

        if (${ITK_FOUND})
                target_link_libraries(meximutil PUBLIC 
                        ${ITK_LIBRARIES})
                target_include_directories(meximutil PUBLIC 
                        ${ITK_INCLUDE_DIRS})
        else ()
                target_link_libraries (meximutil PUBLIC 
                        ${NIFTI_NIFTIIO_LIBRARY})
                target_include_directories (meximutil PUBLIC 
                        ${NIFTI_INCLUDE_DIR})
        endif ()
endif ()

# Configure the installation
install (TARGETS imutil EXPORT SIFT3D-targets 
	 RUNTIME DESTINATION ${INSTALL_BIN_DIR} 
	 LIBRARY DESTINATION ${INSTALL_LIB_DIR} 
	 ARCHIVE DESTINATION ${INSTALL_LIB_DIR})

# Add the code snippets
add_subdirectory (templates)
