/* -----------------------------------------------------------------------------
 * imutil.h
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Public header for imutil.c
 * -----------------------------------------------------------------------------
 */

#include <sys/types.h>
#include "types.h"

#ifndef UTIL_H
#define UTIL_H

/* Vendor-specific info */
#define PLATFORM_NAME_NVIDIA "NVIDIA CUDA"

/* Externally-visible routines */
void clFinish_all();

void check_cl_error(int err, const char *msg);

int init_cl(CL_data *user_cl_data, const char *platform_name, 
			cl_device_type device_type,	cl_mem_flags mem_flags, 
			cl_image_format image_format);

int convert_Mat_rm(const Mat_rm *const in, Mat_rm *const out, 
        const data_type type);

int init_Mat_rm(Mat_rm *mat, int num_rows, int num_cols,
				data_type type, int set_zero);

int init_Mat_rm_p(Mat_rm *mat, const void *p, int num_rows, 
				  int num_cols,	data_type type, int set_zero);

void sprint_type_Mat_rm(const Mat_rm *const mat, char *const str);

int concat_h_Mat_rm(const Mat_rm *const left, const Mat_rm *const right,
        Mat_rm *const dst);

int set_Mat_rm_zero(Mat_rm *mat);

int copy_Mat_rm(const Mat_rm *const src, Mat_rm *const dst);

int print_Mat_rm(Mat_rm *mat);

int resize_Mat_rm(Mat_rm *mat); 

int eigen_Mat_rm(Mat_rm *A, Mat_rm *Q, Mat_rm *L);

int solve_Mat_rm(Mat_rm *A, Mat_rm *B, double limit, Mat_rm *X);

int solve_Mat_rm_ls(Mat_rm *A, Mat_rm *B, Mat_rm *X);

int transpose_Mat_rm(Mat_rm *src, Mat_rm *dst);

int det_symm_Mat_rm(Mat_rm *mat, void *det);

int zero_Mat_rm(Mat_rm *mat);
 
void cleanup_Mat_rm(Mat_rm *mat);

int init_tform(void *tform, const tform_type type);

int init_Affine(Affine *affine, int dim);

int copy_tform(const void *const src, void *const dst);

int Affine_set_mat(const Mat_rm *const mat, Affine *const affine);

void apply_tform_xyz(void *const tform, const double x_in, const double y_in, 
        const double z_in, double *const x_out, double *const y_out, 
        double *const z_out);

int apply_tform_Mat_rm(void *const tform, const Mat_rm *const mat_in, 
        Mat_rm *const mat_out);

tform_type tform_get_type(const void *const tform);

size_t tform_get_size(const void *const tform);

size_t tform_type_get_size(const tform_type type);

void cleanup_tform(void *const tform);

int write_tform(const char *path, const void *const tform);

int mul_Mat_rm(const Mat_rm *const mat_in1, const Mat_rm *const mat_in2, 
        Mat_rm *const mat_out);

int init_Sep_FIR_filter(Sep_FIR_filter *f, int dim, int half_width, int width, 
						float *kernel, int symmetric);

void cleanup_Sep_FIR_filter(Sep_FIR_filter *f);

int apply_Sep_FIR_filter(const Image *const src, Image *const dst, 
        Sep_FIR_filter *const f);

int draw_grid(Image *grid, int nx, int ny, int nz, int spacing, 
					   int line_width);

int draw_points(const Mat_rm *const in, const int *const dims, int radius, 
                Image *const out);

int draw_lines(const Mat_rm *const points1, const Mat_rm *const points2, 
	       const int *const dims, Image *const out);

int read_nii(const char *path, Image *im);

int write_nii(const char *path, Image *im);

int write_Mat_rm(const char *path, const Mat_rm *const mat);

int init_im_first_time(Image *im, const int nx, const int ny, const int nz,
                        const int nc);

int im_load_cl(Image *im, int blocking);

int im_copy_dims(const Image *const src, Image *dst);

int im_copy_data(const Image *const src, Image *const dst);

void im_free(Image *im);

int im_downsample_2x(Image *src, Image *dst);

int im_downsample_2x_cl(Image *src, Image *dst);

int im_read_back(Image *im, int blocking);

int im_set_kernel_arg(cl_kernel kernel, int n, Image *im);

int im_transpose(const Image *const src, const int dim1, const int dim2, 
		 Image *const dst);

int im_upsample_2x(Image *src, Image *dst);

int im_pad(Image *im, Image *pad);

void im_default_stride(Image *im);

int im_resize(Image *im);

int im_concat(const Image *const src1, const Image *const src2, const int dim, 
	      Image *const dst);

void im_scale(Image *im);

int im_subtract(Image *src1, Image *src2, Image *dst);

void im_zero(Image *im);

void im_Hessian(Image *im, int x, int y, int z, Mat_rm *H);

int im_inv_transform(void *const tform, const Image *const in, 
        Image *const out, const interp_type interp);

void init_im(Image *im);

int init_Gauss_filter(Gauss_filter *const gauss, const double sigma, 
                      const int dim);

int init_Gauss_incremental_filter(Gauss_filter *const gauss, 
                const double s_cur, const double s_next, const int dim);

int init_Sep_FIR_filter(Sep_FIR_filter *f, int dim, int half_width, int width, 
						float *kernel, int symmetric);

void cleanup_Sep_FIR_filter(Sep_FIR_filter *f);

void cleanup_Gauss_filter(Gauss_filter *gauss);

void init_GSS_filters(GSS_filters *const gss);

int make_gss(GSS_filters *const gss, const Pyramid *const pyr);

void cleanup_GSS_filters(GSS_filters *const gss);

void init_Pyramid(Pyramid *const pyr);

int copy_Pyramid(const Pyramid *const src, Pyramid *const dst);

int resize_Pyramid(Pyramid *pyr, Image *im);

void cleanup_Pyramid(Pyramid *const pyr);

void init_Slab(Slab *const slab);

void cleanup_Slab(Slab *const slab);

int resize_Slab(Slab *slab, int num, size_t size);

int write_pyramid(const char *path, Pyramid *pyr);

void err_exit(const char *str);

int mkpath(const char *path, mode_t mode);

int init_Ransac(Ransac *ran); 
					  
int set_min_inliers_Ransac(Ransac *ran, double min_inliers);

int set_err_thresh_Ransac(Ransac *ran, double err_thresh);

void set_num_iter_Ransac(Ransac *ran, unsigned int num_iter);

int find_tform_ransac(Ransac* ran, Mat_rm* src, Mat_rm* ref, const int dim,
		      void *const tform);

int parse_gnu(const int argc, char *const *argv);

void print_bug_msg();

#endif
