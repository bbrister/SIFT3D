################################################################################
# Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the Matlab wrapper toolbox.
################################################################################

# There is an issue with MEX compilation for Matlab 2018b that requires a newer 
# CMake.
cmake_minimum_required (VERSION 3.15)

# Enumerate the Matlab scripts for testing and installation
set (INSTALL_SCRIPTS setupSift3D.m imRead3D.m imWrite3D.m detectSift3D.m
        extractSift3D.m orientation3D.m keypoint3D.m registerSift3D.m 
        matchSift3D.m checkUnits3D.m Sift3DParser.m)
set (TEST_SCRIPTS Sift3DTest.m)

# Build the mex utility library
add_library (mexutil SHARED mexutil.c)
target_include_directories (mexutil PUBLIC ${Matlab_INCLUDE_DIRS})
target_include_directories (mexutil PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)
target_link_libraries (mexutil PUBLIC mexreg mexsift3D meximutil
        ${Matlab_LIBRARIES})

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

matlab_add_mex (NAME mexOrientation3D
        SRC mexOrientation3D.c
        LINK_TO mexutil ${Matlab_LIBRARIES}
)

matlab_add_mex (NAME mexExtractSift3D
        SRC mexExtractSift3D.c
        LINK_TO mexutil ${Matlab_LIBRARIES}
)

matlab_add_mex (NAME mexRegisterSift3D
        SRC mexRegisterSift3D.c
        LINK_TO mexutil ${Matlab_LIBRARIES}
)

matlab_add_mex (NAME mexMatchSift3D
        SRC mexMatchSift3D.c
        LINK_TO mexutil ${Matlab_LIBRARIES}
)

# Enumerate the binary targets
set (BINARIES mexutil mexImRead3D mexImWrite3D mexDetectSift3D 
        mexOrientation3D mexExtractSift3D mexRegisterSift3D mexMatchSift3D)

# Configure the build destination
set_target_properties (${BINARIES}
        PROPERTIES 
        ARCHIVE_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
        LIBRARY_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
        RUNTIME_OUTPUT_DIRECTORY ${BUILD_TOOLBOX_DIR}
)

# Copy the scripts and test suite to the build tree
foreach (FILENAME IN LISTS INSTALL_SCRIPTS TEST_SCRIPTS)
        configure_file (${FILENAME} "${BUILD_TOOLBOX_DIR}/${FILENAME}" COPYONLY)
endforeach ()

# Install the matlab scripts
install (FILES ${INSTALL_SCRIPTS} 
        DESTINATION ${INSTALL_TOOLBOX_DIR}
)

# Install the libraries and mex files
install (TARGETS ${BINARIES}
        RUNTIME DESTINATION ${INSTALL_TOOLBOX_DIR}
        LIBRARY DESTINATION ${INSTALL_TOOLBOX_DIR}
        ARCHIVE DESTINATION ${INSTALL_TOOLBOX_DIR}
)
