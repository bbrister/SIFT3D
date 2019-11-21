################################################################################
# Copyright (c) 2019 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the Python package
################################################################################

# Create an installation subdirectory for the python module
set (INSTALL_PYTHON_MODULE_DIR "${INSTALL_PYTHON_DIR}/pySift3D")

# Set the install location for Python scripts
set (INSTALL_SCRIPTS __init__.py imutil.py)
install (FILES ${INSTALL_SCRIPTS} 
        DESTINATION ${INSTALL_PYTHON_MODULE_DIR}
)
