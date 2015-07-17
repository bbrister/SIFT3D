################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for 3D SIFT keypoint detection, description, and matching.
################################################################################

add_library (sift3D SHARED sift.c)
target_link_libraries (sift3D PUBLIC imutil ${M_LIBRARY})
target_include_directories (sift3D PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library (mexsift3D SHARED sift.c)
target_link_libraries (mexsift3D PUBLIC meximutil ${M_LIBRARY})
target_include_directories (mexsift3D PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties (mexsift3D
        PROPERTIES 
        ARCHIVE_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        LIBRARY_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        RUNTIME_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
)

install (TARGETS sift3D
	 RUNTIME DESTINATION bin 
	 LIBRARY DESTINATION lib/sift3d
	 ARCHIVE DESTINATION lib/sift3d)
