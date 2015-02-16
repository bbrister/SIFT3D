######################################################################
# Build file for 3D SIFT keypoint detection, description, and matching
######################################################################

add_library (sift3D SHARED sift.c)
target_link_libraries (sift3D PUBLIC imutil)
target_include_directories (sift3D PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable (getSift3DFeatures driver.c)
target_link_libraries (getSift3DFeatures PUBLIC imutil sift3D)

install (TARGETS sift3D getSift3DFeatures 
	 RUNTIME DESTINATION bin 
	 LIBRARY DESTINATION lib
	 ARCHIVE DESTINATION lib)
