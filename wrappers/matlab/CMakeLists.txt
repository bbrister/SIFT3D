################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the Matlab wrapper toolbox.
################################################################################

# Build the mex utility library
add_library (mexutil SHARED mexutil.c)
target_include_directories (mexutil PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)
target_link_libraries (mexutil PUBLIC mexsift3D meximutil)

# Build the mex files
matlab_add_mex (NAME mexDetectSift3D 
        SRC mexDetectSift3D.c 
        LINK_TO mexutil meximutil -lmwlapack -lmwblas
)

# Copy the Matlab files
file (COPY setupSift3D.m detectSift3D.m 
        DESTINATION ${TOOLBOX_PATH})

# Send all files to the toolbox subdirectory 
set_target_properties (mexDetectSift3D
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        LIBRARY_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        RUNTIME_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
)
