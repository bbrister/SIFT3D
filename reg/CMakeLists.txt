######################################################################
# Build file for feature-based registration 
######################################################################

add_executable (siftreg siftreg.c)
target_link_libraries (siftreg PUBLIC imutil sift3D)

if (UNIX)
    target_link_libraries (siftreg PUBLIC m)
endif ()

install (TARGETS siftreg
	 RUNTIME DESTINATION bin 
	 LIBRARY DESTINATION lib/sift3d
	 ARCHIVE DESTINATION lib/sift3d)
