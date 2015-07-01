################################################################################
# Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
################################################################################
# Build file for the C++/ITK wrappers.
################################################################################

# Skip the build if ITK is not installed
if (ITK_FOUND)
	install (FILES itkwrap.h itkwrap.inl DESTINATION include/sift3d) 
endif ()
