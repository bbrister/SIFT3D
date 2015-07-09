################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for example programs.
################################################################################

# The data directory
set (DATA_DIR "data")

# The name of the examples subdirectory, where compiled examples are found 
set (EXAMPLES_DIR "examples")
set (EXAMPLES_PATH "${CMAKE_BINARY_DIR}/${EXAMPLES_DIR}")

# Copy the example data
configure_file("${DATA_DIR}/1.nii.gz" "${EXAMPLES_PATH}/1.nii.gz" COPYONLY)
configure_file("${DATA_DIR}/2.nii.gz" "${EXAMPLES_PATH}/2.nii.gz" COPYONLY)

# The example programs
add_executable (featuresC featuresC.c)
target_link_libraries (featuresC PUBLIC sift3D imutil)

add_executable (registerC registerC.c)
target_link_libraries (registerC PUBLIC reg sift3D imutil)

# Send all files to the examples subdirectory 
set_target_properties(featuresC registerC
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${EXAMPLES_PATH}
        LIBRARY_OUTPUT_DIRECTORY ${EXAMPLES_PATH}
        RUNTIME_OUTPUT_DIRECTORY ${EXAMPLES_PATH}
)
