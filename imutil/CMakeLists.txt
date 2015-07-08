################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the image processing utility library.
################################################################################

add_library (imutil SHARED imutil.c)
target_include_directories (imutil PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (imutil PUBLIC ${LAPACK_LIBRARIES} ${Z_LIBRARY} ${M_LIBRARY})
install (FILES types.h macros.h imutil.h kernels.cl DESTINATION include/sift3d)

if (ITK_FOUND)
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

install (TARGETS imutil 
	 RUNTIME DESTINATION bin 
	 LIBRARY DESTINATION lib/sift3d
	 ARCHIVE DESTINATION lib/sift3d)

add_subdirectory (templates)
