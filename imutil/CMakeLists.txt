################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the image processing utility library.
################################################################################

# Find LAPACK
find_package (LAPACK REQUIRED)

# Find DCMTK
find_package (DCMTK REQUIRED)

# Compile imutil
add_library (imutil SHARED imutil.c dicom.cpp)
target_include_directories (imutil PUBLIC 
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)
target_link_libraries (imutil PUBLIC ${LAPACK_LIBRARIES} ${DCMTK_LIBRARIES})
install (FILES types.h macros.h imutil.h kernels.cl dicom.h
        DESTINATION ${INSTALL_INCLUDE_DIR})

# Try to find ITK
find_package (ITK 4 QUIET)

if (${ITK_FOUND})
	message(STATUS "Linking to ITK")

	target_link_libraries(imutil PUBLIC ${ITK_LIBRARIES})
	target_include_directories(imutil PUBLIC ${ITK_INCLUDE_DIRS})

# Otherwise link directly to libniftiio
else ()
	message(STATUS "ITK not found. Linking to niftiio directly")
	find_package (NIFTI REQUIRED)

	target_link_libraries (imutil PUBLIC ${NIFTI_NIFTIIO_LIBRARY})
	target_include_directories (imutil PUBLIC ${NIFTI_INCLUDE_DIR})

endif ()

# Link to system libraries
target_link_libraries(imutil PUBLIC ${Z_LIBRARY} ${M_LIBRARY})

# If Matlab was found, compile a copy for use with Matlab libraries
if (${Matlab_FOUND})

        add_library (meximutil SHARED imutil.c dicom.cpp)
        target_compile_definitions (meximutil PUBLIC "-DSIFT3D_MEX")
        set_target_properties (meximutil 
                PROPERTIES 
                ARCHIVE_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
                LIBRARY_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
                RUNTIME_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        )
        target_include_directories (meximutil PUBLIC 
                ${Matlab_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
        target_link_libraries (meximutil PUBLIC ${Z_LIBRARY} ${M_LIBRARY} 
                ${DCMTK_LIBRARIES})

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
