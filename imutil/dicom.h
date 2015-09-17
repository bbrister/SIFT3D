/* -----------------------------------------------------------------------------
 * dicom.h 
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Internal header file for the DCMTK wrapper.
 * -----------------------------------------------------------------------------
 */

#ifndef _DICOM_H
#define _DICOM_H

#ifdef __cplusplus
extern "C" {
#endif


/* Dicom file extension */
const char ext_dcm[] = "dcm";

/* Internal struct to hold limited Dicom metadata */
typedef struct _Dcm_meta {
        const char *patient_name; // Patient name
        const char *patient_id; // Patient ID
        const char *study_uid; // Study Instance UID
        const char *series_uid; // Series UID
        const char *series_descrip;
        const char *instance_uid; // SOP Instance UID
        int instance_num; // Instance number
} Dcm_meta;

int read_dcm(const char *path, Image *const im);

int read_dcm_dir(const char *path, Image *const im);

int write_dcm(const char *path, const Image *const im, 
        const Dcm_meta *const meta);

int write_dcm_dir(const char *path, const Image *const im, 
        const Dcm_meta *const meta);

#ifdef __cplusplus
}
#endif

#endif
