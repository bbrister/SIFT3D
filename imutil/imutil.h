/* -----------------------------------------------------------------------------
 * imutil.h
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2017 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Public header for imutil.c
 * -----------------------------------------------------------------------------
 */

#include "imtypes.h"

#ifndef _IMUTIL_H
#define _IMUTIL_H

#ifdef __cplusplus
extern "C" {
#endif

/* Extra return codes for this module */
#define SIFT3D_FILE_DOES_NOT_EXIST 1 /* The file does not exist */
#define SIFT3D_UNSUPPORTED_FILE_TYPE 2 /* The file type is not supported */
#define SIFT3D_WRAPPER_NOT_COMPILED 3 /* The file type is supported, but the 
                                     * wrapper library was not compiled. */
#define SIFT3D_UNEVEN_SPACING 4 /* The image slices are not evenly spaced. */
#define SIFT3D_INCONSISTENT_AXES 5 /* The image slices have inconsistent 
                                    * axes. */
#define SIFT3D_DUPLICATE_SLICES 6 /* Multiple slices in the same location. */

/* Vendor-specific info */
#define PLATFORM_NAME_NVIDIA "NVIDIA CUDA"

/* Parameters */
const extern double SIFT3D_err_thresh_default;
const extern int SIFT3D_num_iter_default;

/* Externally-visible routines */
void *SIFT3D_safe_realloc(void *ptr, size_t size);

void clFinish_all();

void check_cl_error(int err, const char *msg);

int init_cl(CL_data *user_cl_data, const char *platform_name, 
			cl_device_type device_type,	cl_mem_flags mem_flags, 
			cl_image_format image_format);

void init_Mesh(Mesh * const mesh);

void cleanup_Mesh(Mesh * const mesh);

int convert_Mat_rm(const Mat_rm *const in, Mat_rm *const out, 
        const Mat_rm_type type);

int init_Mat_rm(Mat_rm *const mat, const int num_rows, const int num_cols,
                const Mat_rm_type type, const int set_zero);

int init_Mat_rm_p(Mat_rm *const mat, const void *const p, const int num_rows, 
                  const int num_cols, const Mat_rm_type type, 
                  const int set_zero);

void sprint_type_Mat_rm(const Mat_rm *const mat, char *const str);

int concat_Mat_rm(const Mat_rm * const src1, const Mat_rm * const src2,
		    Mat_rm * const dst, const int dim);

int set_Mat_rm_zero(Mat_rm *mat);

int copy_Mat_rm(const Mat_rm *const src, Mat_rm *const dst);

int print_Mat_rm(const Mat_rm *const mat);

int resize_Mat_rm(Mat_rm *const mat); 

int eigen_Mat_rm(Mat_rm *A, Mat_rm *Q, Mat_rm *L);

int solve_Mat_rm(const Mat_rm *const A, const Mat_rm *const B, 
        const double limit, Mat_rm *const X);

int solve_Mat_rm_ls(const Mat_rm *const A, const Mat_rm *const B, 
        Mat_rm *const X);

int transpose_Mat_rm(const Mat_rm *const src, Mat_rm *const dst);

int det_symm_Mat_rm(Mat_rm *mat, void *det);

int zero_Mat_rm(Mat_rm *const mat);
 
int identity_Mat_rm(const int n, Mat_rm *const mat);

void cleanup_Mat_rm(Mat_rm *mat);

int init_tform(void *const tform, const tform_type type);

int init_Affine(Affine *const affine, const int dim);

int copy_tform(const void *const src, void *const dst);

int Affine_set_mat(const Mat_rm *const mat, Affine *const affine);

void apply_tform_xyz(const void *const tform, const double x_in, 
                     const double y_in, const double z_in, double *const x_out,
		     double *const y_out, double *const z_out);

int apply_tform_Mat_rm(const void *const tform, const Mat_rm *const mat_in, 
        Mat_rm *const mat_out);

tform_type tform_get_type(const void *const tform);

size_t tform_get_size(const void *const tform);

size_t tform_type_get_size(const tform_type type);

void cleanup_tform(void *const tform);

int write_tform(const char *path, const void *const tform);

int mul_Mat_rm(const Mat_rm *const mat_in1, const Mat_rm *const mat_in2, 
        Mat_rm *const mat_out);

int draw_grid(Image *grid, int nx, int ny, int nz, int spacing, 
					   int line_width);

int draw_points(const Mat_rm *const in, const int *const dims, int radius, 
                Image *const out);

int draw_lines(const Mat_rm *const points1, const Mat_rm *const points2, 
	       const int *const dims, Image *const out);

im_format im_get_format(const char *path);

int im_read(const char *path, Image *const im);

int im_write(const char *path, const Image *const im);

char *im_get_parent_dir(const char *path);

int write_Mat_rm(const char *path, const Mat_rm *const mat);

int init_im_with_dims(Image *const im, const int nx, const int ny, const int nz,
                        const int nc);

int im_load_cl(Image *im, int blocking);

int im_copy_dims(const Image *const src, Image *dst);

int im_copy_data(const Image *const src, Image *const dst);

void im_free(Image *im);

int im_channel(const Image * const src, Image * const dst,
	       const unsigned int chan);

int im_downsample_2x(const Image *const src, Image *const dst);

int im_downsample_2x_cl(Image *src, Image *dst);

int im_read_back(Image *im, int blocking);

int im_set_kernel_arg(cl_kernel kernel, int n, Image *im);

int im_permute(const Image *const src, const int dim1, const int dim2, 
		 Image *const dst);

int im_upsample_2x(const Image *const src, Image *const dst);

int im_pad(const Image *const im, Image *const pad);

void im_default_stride(Image *const im);

int im_resize(Image *const im);

int im_concat(const Image *const src1, const Image *const src2, const int dim, 
	      Image *const dst);

float im_max_abs(const Image *const im);

void im_scale(const Image *const im);

int im_subtract(Image *src1, Image *src2, Image *dst);

void im_zero(Image *im);

void im_Hessian(Image *im, int x, int y, int z, Mat_rm *H);

int im_inv_transform(const void *const tform, const Image * const src,
		     const interp_type interp, const int resize, 
                     Image *const dst);

int im_resample(const Image *const src, const double *const units, 
	const interp_type interp, Image *const dst);

void init_im(Image *const im);

int init_Gauss_filter(Gauss_filter *const gauss, const double sigma, 
                      const int dim);

int init_Gauss_incremental_filter(Gauss_filter *const gauss, 
                const double s_cur, const double s_next, const int dim); 

int init_Sep_FIR_filter(Sep_FIR_filter *const f, const int dim, const int width,
			const float *const kernel, const int symmetric);

int apply_Sep_FIR_filter(const Image *const src, Image *const dst, 
        Sep_FIR_filter *const f, const double unit);

void cleanup_Sep_FIR_filter(Sep_FIR_filter *const f);

void cleanup_Gauss_filter(Gauss_filter *gauss);

void init_GSS_filters(GSS_filters *const gss);

int make_gss(GSS_filters *const gss, const Pyramid *const pyr);

void cleanup_GSS_filters(GSS_filters *const gss);

void init_Pyramid(Pyramid *const pyr);

int copy_Pyramid(const Pyramid *const src, Pyramid *const dst);

int resize_Pyramid(const Image *const im, const int first_level, 
        const unsigned int num_kp_levels, const unsigned int num_levels,
        const int first_octave, const unsigned int num_octaves, 
        Pyramid *const pyr);

int set_scales_Pyramid(const double sigma0, const double sigma_n, 
        Pyramid *const pyr);

void cleanup_Pyramid(Pyramid *const pyr);

void init_Slab(Slab *const slab);

void cleanup_Slab(Slab *const slab);

int resize_Slab(Slab *slab, int num, size_t size);

int write_pyramid(const char *path, Pyramid *pyr);

void err_exit(const char *str);

void init_Ransac(Ransac *const ran);
					  
int set_err_thresh_Ransac(Ransac *const ran, double err_thresh);

int set_num_iter_Ransac(Ransac *const ran, int num_iter);

int copy_Ransac(const Ransac *const src, Ransac *const dst);

int find_tform_ransac(const Ransac *const ran, const Mat_rm *const src, 
        const Mat_rm *const ref, void *const tform);

int parse_gnu(const int argc, char *const *argv);

void print_bug_msg();

#ifdef __cplusplus
}
#endif

#endif
