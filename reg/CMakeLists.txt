################################################################################
# Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for image registration.
################################################################################

add_library (reg SHARED reg.c)
target_link_libraries (reg PUBLIC sift3D imutil)
target_link_libraries (reg PRIVATE ${M_LIBRARY})
target_include_directories (reg PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)
install (FILES reg.h DESTINATION ${INSTALL_INCLUDE_DIR})

install (TARGETS reg EXPORT SIFT3D-targets
	 RUNTIME DESTINATION ${INSTALL_BIN_DIR} 
	 LIBRARY DESTINATION ${INSTALL_LIB_DIR} 
	 ARCHIVE DESTINATION ${INSTALL_LIB_DIR})

# If Matlab was found, compile a copy for use with matlab wrappers
if (BUILD_Matlab)

        add_library (mexreg SHARED reg.c)
        target_link_libraries (mexreg PUBLIC mexsift3D meximutil)
        target_link_libraries (mexreg PRIVATE ${M_LIBRARY})
        target_include_directories (mexreg PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
        )
        set_target_properties (mexreg
                PROPERTIES 
                ARCHIVE_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
                LIBRARY_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
                RUNTIME_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
        )

        install (TARGETS mexreg
                RUNTIME DESTINATION ${INSTALL_TOOLBOX_DIR}
                LIBRARY DESTINATION ${INSTALL_TOOLBOX_DIR}
                ARCHIVE DESTINATION ${INSTALL_TOOLBOX_DIR}
        )
endif ()
