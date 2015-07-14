################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the Matlab wrapper toolbox.
################################################################################

# The name and location of the Matlab toolbox
set (TOOLBOX_NAME "matlab")
set (TOOLBOX_PATH "${WRAPPERS_PATH}/${TOOLBOX_NAME}")

# The mex files are shared libraries
matlab_add_mex (NAME mexDetectSift3D 
        SRC mexDetectSift3D.c 
        LINK_TO sift3D imutil)

# Copy the Matlab files
configure_file (detectSift3D.m ${TOOLBOX_PATH} COPYONLY)

# Send all files to the toolbox subdirectory 
set_target_properties (mexDetectSift3D
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        LIBRARY_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        RUNTIME_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
)
