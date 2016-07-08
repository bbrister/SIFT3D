################################################################################
# Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the command line interfaces.
################################################################################

add_executable (denseSift3D denseSift3D.c)
target_link_libraries(denseSift3D PUBLIC sift3D imutil)
target_link_libraries(denseSift3D PRIVATE ${M_LIBRARY})

add_executable(kpSift3D kpSift3D.c)
target_link_libraries(kpSift3D PUBLIC sift3D imutil)

add_executable(regSift3D regSift3D.c)
target_link_libraries(regSift3D PUBLIC reg sift3D imutil)

install (TARGETS denseSift3D kpSift3D regSift3D
	 RUNTIME DESTINATION ${INSTALL_BIN_DIR} 
	 LIBRARY DESTINATION ${INSTALL_LIB_DIR} 
	 ARCHIVE DESTINATION ${INSTALL_LIB_DIR})
