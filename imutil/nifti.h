/* -----------------------------------------------------------------------------
 * nifti.h 
 * -----------------------------------------------------------------------------
 * Copyright (c) 2017 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Internal header file for the nifticlib wrapper.
 * -----------------------------------------------------------------------------
 */

#ifndef _NIFTI_H
#define _NIFTI_H

int read_nii(const char *path, Image *const im);

int write_nii(const char *path, const Image *const im);

#endif
