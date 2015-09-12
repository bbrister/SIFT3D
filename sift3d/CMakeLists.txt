################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for 3D SIFT keypoint detection, description, and matching.
################################################################################

add_library (sift3D SHARED sift.c)
target_link_libraries (sift3D PUBLIC imutil ${M_LIBRARY})
target_include_directories (sift3D PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)
install (FILES sift.h DESTINATION ${INSTALL_INCLUDE_DIR})

# If Matlab was found, compile a copy for use with matlab wrappers
if (${Matlab_FOUND})

        add_library (mexsift3D SHARED sift.c)
        target_include_directories (mexsift3D PUBLIC 
                ${CMAKE_CURRENT_SOURCE_DIR})
        target_link_libraries (mexsift3D PUBLIC meximutil ${M_LIBRARY})
        set_target_properties (mexsift3D
                PROPERTIES 
                ARCHIVE_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
                LIBRARY_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
                RUNTIME_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
        )
        install (TARGETS mexsift3D
                RUNTIME DESTINATION ${INSTALL_TOOLBOX_DIR}
                LIBRARY DESTINATION ${INSTALL_TOOLBOX_DIR}
                ARCHIVE DESTINATION ${INSTALL_TOOLBOX_DIR}
        )
endif ()

install (TARGETS sift3D EXPORT SIFT3D-targets
	 RUNTIME DESTINATION ${INSTALL_BIN_DIR} 
	 LIBRARY DESTINATION ${INSTALL_LIB_DIR} 
	 ARCHIVE DESTINATION ${INSTALL_LIB_DIR})
