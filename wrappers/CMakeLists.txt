################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the wrappers for other programming languages.
################################################################################

# Subdirectories for each language
if (${Matlab_FOUND})
        add_subdirectory (matlab)
endif ()
