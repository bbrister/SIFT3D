################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for image registration.
################################################################################

add_library (reg SHARED reg.c)
target_link_libraries (reg PUBLIC sift3D imutil)
target_include_directories (reg PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)
install (FILES reg.h DESTINATION ${INSTALL_INCLUDE_DIR})

install (TARGETS reg EXPORT SIFT3D-targets
	 RUNTIME DESTINATION ${INSTALL_BIN_DIR} 
	 LIBRARY DESTINATION ${INSTALL_LIB_DIR} 
	 ARCHIVE DESTINATION ${INSTALL_LIB_DIR})
