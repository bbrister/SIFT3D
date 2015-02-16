######################################################################
# Build file for the image processing utility library
######################################################################

add_library (imutil SHARED imutil.c)
target_include_directories (imutil PUBLIC $(CMAKE_CURRENT_SOURCE_DIR))
install (FILES types.h macros.h imutil.h kernels.cl DESTINATION include/sift3d)

# Add the dependencies
# XXX replace with ITK
target_link_libraries (imutil PUBLIC lapack niftiio)
target_include_directories (imutil PUBLIC /usr/include/nifti)

install (TARGETS imutil 
	 RUNTIME DESTINATION bin 
	 LIBRARY DESTINATION lib
	 ARCHIVE DESTINATION lib)

add_subdirectory (templates)
