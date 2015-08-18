/* -----------------------------------------------------------------------------
 * dicom.hpp 
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Header file for the DCMTK wrapper.
 * -----------------------------------------------------------------------------
 */

#ifndef _DICOM_HPP
#define _DICOM_HPP

#ifdef __cplusplus
extern "C" {
#endif

int read_dcm(const char *path, Image *const im);

int read_dcm_dir(const char *path, Image *const im);

int write_dcm(const char *path, const Image *const im);

int write_dcm_dir(const char *path, const Image *const im);

#ifdef __cplusplus
}
#endif

#endif
