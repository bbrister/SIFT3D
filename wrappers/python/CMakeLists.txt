################################################################################
# Copyright (c) 2019 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the Python wrapper module.
################################################################################

# Compile the C compatibility layer
add_library (sift3Dpython SHARED sift3dpython.c)
target_link_libraries (sift3Dpython PRIVATE imutil)

# Install the python setup script
set (INSTALL_SCRIPTS setup.py)
install (FILES ${INSTALL_SCRIPTS}
        DESTINATION ${INSTALL_PYTHON_DIR}
)

# Install the C library
install (
        TARGETS sift3Dpython 
        EXPORT SIFT3D-targets
        RUNTIME DESTINATION ${INSTALL_BIN_DIR} 
        LIBRARY DESTINATION ${INSTALL_LIB_DIR} 
        ARCHIVE DESTINATION ${INSTALL_LIB_DIR}
)

# Include the python scipts subdirectory
add_subdirectory (pySift3D)
