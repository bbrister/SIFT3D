################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the Matlab wrapper toolbox.
################################################################################

# Build the mex utility library
add_library (mexutil SHARED mexutil.c)
target_include_directories (mexutil PUBLIC ${Matlab_INCLUDE_DIRS})
target_include_directories (mexutil PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)
target_link_libraries (mexutil PUBLIC mexsift3D meximutil 
        ${Matlab_MWLAPACK_LIBRARY} ${Matlab_MWBLAS_LIBRARY} ${Matlab_LIBRARIES})

# Build the mex files
matlab_add_mex (NAME mexImRead3D
        SRC mexImRead3D.c
        LINK_TO mexutil ${Matlab_LIBRARIES}
)

matlab_add_mex (NAME mexImWrite3D
        SRC mexImWrite3D.c
        LINK_TO mexutil ${Matlab_LIBRARIES}
)

matlab_add_mex (NAME mexDetectSift3D 
        SRC mexDetectSift3D.c 
        LINK_TO mexutil ${Matlab_LIBRARIES}
)

matlab_add_mex (NAME mexExtractSift3D
        SRC mexExtractSift3D.c
        LINK_TO mexutil ${Matlab_LIBRARIES}
)

# Configure the build destination
set_target_properties (mexutil mexImRead3D mexImWrite3D mexDetectSift3D
        mexExtractSift3D
        PROPERTIES 
        ARCHIVE_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
        LIBRARY_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
        RUNTIME_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
)

# Copy the scripts and test suite to the build tree
file (COPY setupSift3D.m imRead3D.m imWrite3D.m detectSift3D.m extractSift3D.m
        keypoint3D.m checkUnits3D.m Sift3DTest.m
        DESTINATION ${BUILD_TOOLBOX_DIR}
)

# Install the matlab scripts
install (FILES setupSift3D.m imRead3D.m imWrite3D.m detectSift3D.m 
        extractSift3D.m keypoint3D.m checkUnits3D.m
        DESTINATION ${INSTALL_TOOLBOX_DIR}
)

# Install the libraries and mex files
install (TARGETS mexutil mexImRead3D mexImWrite3D mexDetectSift3D 
        mexExtractSift3D
        RUNTIME DESTINATION ${INSTALL_TOOLBOX_DIR}
        LIBRARY DESTINATION ${INSTALL_TOOLBOX_DIR}
        ARCHIVE DESTINATION ${INSTALL_TOOLBOX_DIR}
)
