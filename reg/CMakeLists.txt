################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for image registration.
################################################################################

add_library (reg SHARED reg.c)
target_link_libraries (reg PUBLIC sift3D imutil)
target_include_directories (reg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install (TARGETS reg
	 RUNTIME DESTINATION bin 
	 LIBRARY DESTINATION lib/sift3d
	 ARCHIVE DESTINATION lib/sift3d)
