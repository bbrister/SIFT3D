################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the Matlab wrapper toolbox.
################################################################################

# Get the path of the Matlab libraries
get_filename_component (Matlab_LIBRARIES_PATH ${Matlab_MEX_LIBRARY} DIRECTORY)

# Find the Matlab LAPACK and BLAS libraries
find_library (MWLAPACK_LIBRARY mwlapack PATHS ${Matlab_LIBRARIES_PATH})
find_library (MWBLAS_LIBRARY mwblas PATHS ${Matlab_LIBRARIES_PATH})

# Build the mex utility library
add_library (mexutil SHARED mexutil.c)
target_include_directories (mexutil PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)
target_include_directories(mexutil PUBLIC ${Matlab_INCLUDE_DIR})
target_link_libraries (mexutil PUBLIC mexsift3D meximutil ${MWLAPACK_LIBRARY}
        ${MWBLAS_LIBRARY} ${Matlab_LIBRARIES})

# Build the mex files
matlab_add_mex (NAME mexDetectSift3D 
        SRC mexDetectSift3D.c 
        LINK_TO mexutil ${Matlab_LIBRARIES}
)

matlab_add_mex (NAME mexExtractSift3D
        SRC mexExtractSift3D.c
        LINK_TO mexutil ${Matlab_LIBRARIES}
)

# Copy the Matlab files
file (COPY setupSift3D.m detectSift3D.m extractSift3D.m Sift3DTest.m
        DESTINATION ${TOOLBOX_PATH})

# Send all files to the toolbox subdirectory 
set_target_properties (mexDetectSift3D mexExtractSift3D
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        LIBRARY_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
        RUNTIME_OUTPUT_DIRECTORY ${TOOLBOX_PATH}
)
