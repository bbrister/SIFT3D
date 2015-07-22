################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the image processing utility library.
################################################################################

# Find LAPACK
find_package (LAPACK REQUIRED)

# Compile imutil
add_library (imutil SHARED imutil.c)
target_include_directories (imutil PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (imutil PUBLIC ${LAPACK_LIBRARIES} ${Z_LIBRARY} ${M_LIBRARY})
install (FILES types.h macros.h imutil.h kernels.cl DESTINATION include/sift3d)

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

install (TARGETS imutil 
	 RUNTIME DESTINATION bin 
	 LIBRARY DESTINATION lib/sift3d
	 ARCHIVE DESTINATION lib/sift3d)

add_subdirectory (templates)
