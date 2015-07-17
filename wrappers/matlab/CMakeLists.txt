################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the Matlab wrapper toolbox.
################################################################################

# The mex files are shared libraries
matlab_add_mex (NAME mexDetectSift3D 
        SRC mexDetectSift3D.c 
        LINK_TO mexsift3D meximutil -lmwlapack -lmwblas
)

# Copy the Matlab files
configure_file (detectSift3D.m ${TOOLBOX_PATH} COPYONLY)

# Send all files to the toolbox subdirectory 
set_target_properties (mexDetectSift3D
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        LIBRARY_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        RUNTIME_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
)
