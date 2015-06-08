######################################################################
# Build file for the command line interface
######################################################################

add_executable (denseSift3D denseSift3D.c)
target_link_libraries(denseSift3D PUBLIC sift3D imutil m)

add_executable (diffeoSift3D diffeoSift3D.c)
target_link_libraries(diffeoSift3D PUBLIC reg sift3D imutil m)

install (TARGETS denseSift3D diffeoSift3D
	 RUNTIME DESTINATION bin 
	 LIBRARY DESTINATION lib/sift3d
	 ARCHIVE DESTINATION lib/sift3d)
