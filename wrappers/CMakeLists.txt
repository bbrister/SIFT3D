################################################################################
# Copyright (c) 2015-2019 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the wrappers for other programming languages.
################################################################################

# Subdirectories for each language
if (BUILD_Matlab)
        add_subdirectory (matlab)
endif ()

if (BUILD_Python)
        add_subdirectory (python)
endif ()
