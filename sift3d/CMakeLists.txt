######################################################################
# Build file for 3D SIFT keypoint detection, description, and matching
######################################################################

add_library (sift3D SHARED sift.c)
target_link_libraries (sift3D PUBLIC imutil)
target_include_directories (sift3D PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install (TARGETS sift3D
	 RUNTIME DESTINATION bin 
	 LIBRARY DESTINATION lib/sift3d
	 ARCHIVE DESTINATION lib/sift3d)
