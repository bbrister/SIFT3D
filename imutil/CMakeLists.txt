######################################################################
# Build file for the image processing utility library
######################################################################

add_library (imutil SHARED imutil.c)
target_include_directories (imutil PUBLIC $(CMAKE_CURRENT_SOURCE_DIR))
install (FILES types.h macros.h imutil.h kernels.cl DESTINATION include/sift3d)

# Link to ITK if it is present
find_package(ITK 4)
if (ITK_FOUND)

# Otherwise link directly to liblapack and libniftiio
else ()
	target_link_libraries (imutil PUBLIC lapack niftiio)
	target_include_directories (imutil PUBLIC /usr/include/nifti)
endif ()

install (TARGETS imutil 
	 RUNTIME DESTINATION bin 
	 LIBRARY DESTINATION lib/sift3d
	 ARCHIVE DESTINATION lib/sift3d)

add_subdirectory (templates)
