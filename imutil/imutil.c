/* -----------------------------------------------------------------------------
 * imutil.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2017 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * Miscellaneous utility routines for image processing, linear algebra, and 
 * statistical regression. This library completely defines the Image,
 * Mat_rm, and Ransac types, among others, and stands apart from the other
 * source.
 * -----------------------------------------------------------------------------
 */

#include <assert.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stddef.h>
#include <float.h>
#include <zlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "immacros.h"
#include "imtypes.h"
#include "dicom.h"
#include "nifti.h"
#include "imutil.h"

/* Check for a version number */
#if !defined(SIFT3D_VERSION_NUMBER)
#error imutil.c: Must define the preprocessor macro SIFT3D_VERSION_NUMBER
#endif

/* Stringify a macro name */
#define STR(x) #x

/* Stringify the result of a macro expansion */
#define XSTR(x) STR(x)

/* zlib definitions */
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#include <fcntl.h>
#include <io.h>
#define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#define SET_BINARY_MODE(file)
#endif

/* Implementation parameters */
//#define SIFT3D_USE_OPENCL // Use OpenCL acceleration
#define SIFT3D_RANSAC_REFINE	// Use least-squares refinement in RANSAC

/* Implement strnlen, if it's missing */
#ifndef SIFT3D_HAVE_STRNLEN
size_t strnlen(const char *string, size_t maxlen) {
  const char *end = memchr (string, '\0', maxlen);
  return end ? end - string : maxlen;
}
#endif

/* Implement strndup, if it's missing */
#ifndef SIFT3D_HAVE_STRNDUP
char *strndup(const char *s, size_t n) {
  size_t len = strnlen (s, n);
  char *new = malloc (len + 1);
  if (new == NULL)
    return NULL;
  new[len] = '\0';
  return memcpy (new, s, len);
}
#endif

/* SIFT3D version message */
const char version_msg[] =
    "SIFT3D version " XSTR(SIFT3D_VERSION_NUMBER)  " \n"
    "\n"
    "Source code available at https://github.com/bbrister/SIFT3D\n"
    "\n"
    "Please contact Blaine Rister (blaine@stanford.edu) with questions or "
    "concerns. \n";

/* Bug message */
const char bug_msg[] =
    "SIFT3D has encountered an unexpected error. We would appreciate it \n"
    "if you would report this issue at the following page: \n"
    "       https://github.com/bbrister/SIFT3D/issues \n";

/* Supported file extensions */
extern const char ext_dcm[]; // dicom.h
const char ext_analyze[] = "img";;
const char ext_gz[] = "gz";
const char ext_nii[] = "nii";
const char ext_dir[] = "";

/* Output file permissions */
const mode_t out_mode = 0755;

/* Default parameters */
const double SIFT3D_err_thresh_default = 5.0;
const int SIFT3D_num_iter_default = 500;

/* Declarations for the virtual function implementations */
static int copy_Affine(const void *const src, void *const dst);
static int copy_Tps(const void *const src, void *const dst);
static void apply_Affine_xyz(const void *const affine, const double x_in,
			     const double y_in, const double z_in,
			     double *const x_out, double *const y_out,
			     double *const z_out);
static void apply_Tps_xyz(const void *const tps, const double x_in, 
                          const double y_in, const double z_in, 
                          double *const x_out, double *const y_out, 
                          double *const z_out);
static int apply_Affine_Mat_rm(const void *const affine, 
        const Mat_rm * const mat_in, Mat_rm * const mat_out);
static int apply_Tps_Mat_rm(const void *const tps, const Mat_rm * const mat_in,
			    Mat_rm * const mat_out);
static size_t Affine_get_size(void);
static size_t Tps_get_size(void);
static int write_Affine(const char *path, const void *const tform);
static int write_Tps(const char *path, const void *const tform);
static void cleanup_Affine(void *const affine);
static void cleanup_Tps(void *const tps);
static int mkpath(const char *path, mode_t mode);

/* Virtual function tables */
const Tform_vtable Affine_vtable = {
	copy_Affine,
	apply_Affine_xyz,
	apply_Affine_Mat_rm,
	Affine_get_size,
	write_Affine,
	cleanup_Affine
};

const Tform_vtable Tps_vtable = {
	copy_Tps,
	apply_Tps_xyz,
	apply_Tps_Mat_rm,
	Tps_get_size,
	write_Tps,
	cleanup_Tps
};

/* Internal macros */
#define TFORM_GET_VTABLE(arg) (((Affine *) arg)->tform.vtable)
#define AFFINE_GET_DIM(affine) ((affine)->A.num_rows)

/* Global data */
CL_data cl_data;

/* LAPACK declarations */
#ifdef SIFT3D_MEX
// Set the integer width to Matlab's defined width
#include <uchar.h>
#include "mex.h"
typedef mwSignedIndex fortran_int;
#ifdef _WINDOWS
// Remove underscores from FORTRAN functions
#define dlange_ dlange
#define dgecon_ dgecon
#define dgelss_ dgelss
#define dgetrf_ dgetrf
#define dgetrs_ dgetrs
#define dsyevd_ dsyevd
#endif
#else
typedef int32_t fortran_int;
#endif
extern double dlange_(const char *, const fortran_int *, const fortran_int *,
		      const double *, const fortran_int *, double *);
extern void dgecon_(const char *, const fortran_int *, double *,
		    const fortran_int *, const double *, double *,
		    double *, fortran_int *, fortran_int *);

extern void dgelss_(const fortran_int *, const fortran_int *,
		    const fortran_int *, const double *, const fortran_int *,
		    double *, const fortran_int *, double *, const double *,
		    fortran_int *, double *, const fortran_int *,
		    fortran_int *);

extern void dgetrf_(const fortran_int *, const fortran_int *, double *,
		    const fortran_int *, fortran_int *, fortran_int *);

extern void dgetrs_(const char *, const fortran_int *, const fortran_int *,
		    const double *, const fortran_int *, fortran_int *,
		    double *, const fortran_int *, fortran_int *);

extern void dsyevd_(const char *, const char *, const fortran_int *, double *,
		    const fortran_int *, double *, double *,
		    const fortran_int *, fortran_int *, const fortran_int *,
		    fortran_int *);

/* Internal helper routines */
static char *read_file(const char *path);
static int do_mkdir(const char *path, mode_t mode);
static int cross_mkdir(const char *path, mode_t mode);
static double resample_linear(const Image * const in, const double x,
			      const double y, const double z, const int c);
static double resample_lanczos2(const Image * const in, const double x,
				const double y, const double z, const int c);
static double lanczos(double x, double a);
static int check_cl_image_support(cl_context context, cl_mem_flags mem_flags,
				  cl_image_format image_format,
				  cl_mem_object_type image_type);
static int compile_cl_program_from_source(cl_program * program,
					  cl_context context,
					  cl_device_id * devices,
					  int num_devices, char **src,
					  int num_str);
static int n_choose_k(const int n, const int k, int **ret);
static int make_spline_matrix(Mat_rm * src, Mat_rm * src_in, Mat_rm * sp_src,
			      int K_terms, int *r, int dim);
static int make_affine_matrix(const Mat_rm *const pts_in, const int dim, 
        Mat_rm *const mat_out);
static Mat_rm *extract_ctrl_pts(void *tform, tform_type type);
static Mat_rm *extract_ctrl_pts_Tps(Tps * tps);
static int solve_system(const Mat_rm *const src, const Mat_rm *const ref, 
        void *const tform);
static double tform_err_sq(const void *const tform, const Mat_rm *const src, 
        const Mat_rm *const ref, const int i);
static int ransac(const Mat_rm *const src, const Mat_rm *const ref, 
        const Ransac *const ran, void *tform, int **const cset, int *const len);
static int convolve_sep(const Image * const src,
			Image * const dst, const Sep_FIR_filter * const f,
			const int dim, const double unit);
static int convolve_sep_gen(const Image * const src,
			Image * const dst, const Sep_FIR_filter * const f,
			const int dim, const double unit);
static int convolve_sep_cl(const Image * const src,
			Image * const dst, const Sep_FIR_filter * const f,
			int dim, const double unit);
static int convolve_sep_sym(const Image * const src, Image * const dst,
			    const Sep_FIR_filter * const f, const int dim,
                            const double unit);
static const char *get_file_name(const char *path);
static const char *get_file_ext(const char *name);

/* Unfinished public routines */
int init_Tps(Tps * tps, int dim, int terms);
int resize_Tps(Tps * tps, int num_pts, int dim);

/* As realloc, but frees the underlying pointer and returns NULL on error, or
 * if size is 0 and ptr is non-NULL. */
void *SIFT3D_safe_realloc(void *ptr, size_t size) {

	void *ret;

	// Call realloc and handle failures
	if (size == 0 || (ret = realloc(ptr, size)) == NULL) {
		if (ptr != NULL) {
			free(ptr);
		}
		return NULL;
	}

	return ret;

}

/* Finish all OpenCL command queues. */
void clFinish_all()
{
#ifdef SIFT3D_USE_OPENCL
	int i;
	for (i = 0; i < cl_data.num_devices; i++) {
		clFinish(cl_data.queues[i]);
	}
#endif
}

/* Check the error code and exit on error, unless NDEBUG is
 * defined. If exiting, prints the error type and the cause,
 * given by the msg string. */
void check_cl_error(int err, const char *msg)
{
#ifdef NDEBUG
	return;
#endif
	switch (err) {
	case CL_SUCCESS:
		return;
	default:
		printf("unknown OpenCL error %d \n", err);
	}
	printf("Exiting due to error in: %s \n", msg);
	exit(1);
}

/* Returns SIFT3D_SUCCESS if the specified format is supported for this context,
 * or SIFT3D_FAILURE if it is not. */
SIFT3D_IGNORE_UNUSED
static int check_cl_image_support(cl_context context, cl_mem_flags mem_flags,
				  cl_image_format image_format,
				  cl_mem_object_type image_type)
{
#ifdef SIFT3D_USE_OPENCL
	cl_image_format *image_formats;
	cl_int err;
	cl_uint num_image_formats;
	int i, support;

	err = clGetSupportedImageFormats(context, mem_flags, image_type,
					 0, NULL, &num_image_formats);
	if ((image_formats =
	     malloc(num_image_formats * sizeof(cl_image_format))) == NULL)
		return SIFT3D_FAILURE;
	err |= clGetSupportedImageFormats(context, mem_flags, image_type,
					  num_image_formats, image_formats,
					  NULL);
	check_cl_error(err, "2D image formats");
	support = SIFT3D_FALSE;
	for (i = 0; i < num_image_formats; i++) {
		if (memcmp
		    (image_formats + i, &image_format,
		     sizeof(cl_image_format))) {
			support = SIFT3D_TRUE;
			break;
		}
	}

	return (support == SIFT3D_TRUE) ? SIFT3D_SUCCESS : SIFT3D_FAILURE;
#else
	printf
	    ("check_cl_image_support: This version was not compiled with OpenCL!\n");
	return SIFT3D_FAILURE;
#endif
}

/* Initialize the relevant OpenCL data, using the specified device type,
 * memory flags, and image format. If no such devices are found, or these
 * settings are not supported on the device, returns SIFT3D_FAILURE. Returns SIFT3D_SUCCESS
 * when userData is initialized. 
 *
 * This library saves a copy of user_cl_data for use with future calls. To change 
 * the settings of the library, call init_cl again. */

int init_cl(CL_data * user_cl_data, const char *platform_name,
	    cl_device_type device_type, cl_mem_flags mem_flags,
	    cl_image_format image_format)
{
#ifdef SIFT3D_USE_OPENCL
	Kernels kernels;
	cl_platform_id *platforms;
	cl_device_id *devices;
	char *name, *src;
	cl_context_properties properties[3];
	cl_command_queue *queues;
	cl_program program;
	cl_platform_id platform;
	cl_context context;
	cl_uint num_platforms, num_devices;
	cl_int err;
	cl_bool support;
	size_t size;
	int i;

	// Initialize persistent arrays to NULL
	devices = NULL;
	queues = NULL;
	src = NULL;

	// Query the available platforms and select the specified one
	err = clGetPlatformIDs(0, NULL, &num_platforms);
	if (err != CL_SUCCESS || num_platforms < 1)
		goto init_cl_quit;

	if ((platforms = (cl_platform_id *) malloc(num_platforms *
						   sizeof(cl_platform_id))) ==
	    NULL)
		goto init_cl_quit;
	clGetPlatformIDs(num_platforms, platforms, NULL);
	name = NULL;
	platform = (cl_platform_id) NULL;
	for (i = 0; i < num_platforms; i++) {
		err = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 0,
					NULL, &size);
		name = (char *)realloc(name, size * sizeof(char));
		err = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME,
					size, name, NULL);
		if (!strcmp(name, platform_name)) {
			platform = platforms[i];
		}
	}
	free(platforms);
	free(name);

	if (platform == (cl_platform_id) NULL) {
		printf("init_cl: Failed to find platform %s \n", platform_name);
		goto init_cl_quit;
	}
	// Get the number of devices of the specified type
	devices = NULL;
	err = clGetDeviceIDs(platform, device_type, 0, NULL, &num_devices);
	if (err != CL_SUCCESS && err != CL_DEVICE_NOT_FOUND) {
		check_cl_error(err, "Create context");
	} else if (num_devices > 0) {
		devices =
		    (cl_device_id *) malloc(num_devices * sizeof(cl_device_id));
		err =
		    clGetDeviceIDs(platform, device_type, num_devices, devices,
				   NULL);
		check_cl_error(err, "Get devices");
	}

	if (num_devices <= 0 || devices == NULL) {
		puts("init_cl: No OpenCL devices available \n");
		goto init_cl_quit;
	}
	//TODO: Multiple GPUs on one context does not seem to work. Maybe try multiple
	//      contexts, one per GPU?
	num_devices = 1;

	// Create the context
	properties[0] = CL_CONTEXT_PLATFORM;
	properties[1] = (cl_context_properties) platform;
	properties[2] = 0;
	context =
	    clCreateContext(properties, num_devices, devices, NULL, NULL, &err);
	check_cl_error(err, "Create context \n");

	// Create a command queue on each device 
	if ((queues = malloc(num_devices * sizeof(cl_command_queue))) == NULL)
		goto init_cl_quit;

	for (i = 0; i < num_devices; i++) {
		cl_device_type type;
		err = clGetDeviceInfo(devices[i], CL_DEVICE_TYPE,
				      sizeof(cl_device_type), &type, NULL);
		if (type == device_type) {
			if ((queues[i] = clCreateCommandQueue(context,
							      devices[i], 0,
							      NULL)) == NULL) {
				goto init_cl_quit;
			}
		}
	}

	// Query for support of the desired image and memory format on these devices
	for (i = 0; i < num_devices; i++) {

		// Query for image support
		support = CL_FALSE;
		clGetDeviceInfo(devices[i], CL_DEVICE_IMAGE_SUPPORT,
				sizeof(cl_bool), &support, NULL);
		if (support != CL_TRUE) {
			printf
			    ("init_cl: Images are not supported by device %d \n",
			     i);
			goto init_cl_quit;
		}
		// Query for support of the specified image format in both 2D and 3D images
		check_cl_image_support(context, mem_flags, image_format,
				       CL_MEM_OBJECT_IMAGE2D);
		check_cl_image_support(context, mem_flags, image_format,
				       CL_MEM_OBJECT_IMAGE3D);
	}

	// Load the kernels from a file
	if ((src = read_file(KERNELS_PATH)) == NULL) {
		printf("init_cl: Error reading kernel source file %s",
		       KERNELS_PATH);
		goto init_cl_quit;
	}
	// Compile static programs
	if (compile_cl_program_from_source
	    (&program, context, devices, num_devices, &src, 1))
		goto init_cl_quit;
	kernels.downsample_2x_3d =
	    clCreateKernel(program, "downsample_2x_3d", &err);
	check_cl_error(err, "init_cl: create kernels");
	clReleaseProgram(program);

	// Save data to the user and library copies
	cl_data.platform = platform;
	cl_data.devices = devices;
	cl_data.context = context;
	cl_data.queues = queues;
	cl_data.num_devices = num_devices;
	cl_data.image_format = image_format;
	cl_data.valid = SIFT3D_TRUE;
	cl_data.kernels = kernels;
	*user_cl_data = cl_data;
	return SIFT3D_SUCCESS;

 init_cl_quit:
	if (devices != NULL)
		free(devices);
	if (queues != NULL)
		free(queues);
	if (src != NULL)
		free(src);
	return SIFT3D_FAILURE;
#else
	printf("init_cl: This version was not compiled with OpenCL!\n");
	return SIFT3D_FAILURE;
#endif
}

/* Compile a program from the given source strings, writing the program handle into
 * the specified pointer. */
SIFT3D_IGNORE_UNUSED
static int compile_cl_program_from_source(cl_program * program,
					  cl_context context,
					  cl_device_id * devices,
					  int num_devices, char **src,
					  int num_str)
{
#ifdef SIFT3D_USE_OPENCL
	cl_int err;
	int i;

	err = CL_SUCCESS;
	*program =
	    clCreateProgramWithSource(context, 1, (const char **)src, NULL,
				      &err);
	if (*program == NULL || err != CL_SUCCESS) {
		puts("Error creating program for static kernels \n");
		return SIFT3D_FAILURE;
	}
	if ((err =
	     clBuildProgram(*program, 0, NULL, NULL, NULL,
			    NULL)) != CL_SUCCESS) {
		char log[1 << 15];
		puts("Error: Failed to build program \n");
		for (i = 0; i < num_devices; i++) {
			clGetProgramBuildInfo(*program, devices[i],
					      CL_PROGRAM_BUILD_LOG, sizeof(log),
					      log, NULL);
			printf("\n-------Build log for device %d-------\n %s",
			       i, log);
		}
		return SIFT3D_FAILURE;
	}

	return SIFT3D_SUCCESS;
#else
	printf
	    ("compile_cl_program_from_source: This verison was not compiled with "
	     "OpenCL!\n");
	return SIFT3D_FAILURE;
#endif
}

/* Initialize a triangle mesh for first use. This must be called before mesh
 * can be used in any other functions. */
void init_Mesh(Mesh * const mesh)
{
	mesh->tri = NULL;
	mesh->num = -1;
}

/* Release all memory associated with a triangle mesh. mesh cannot be reused
 * before it is reinitialized. */
void cleanup_Mesh(Mesh * const mesh)
{
	free(mesh->tri);
}

/* Convert a matrix to a different type. in and out may be the same pointer.
 * 
 * This function resizes out.
 * 
 * All matrices must be initialized prior to calling this function. */
int convert_Mat_rm(const Mat_rm * const in, Mat_rm * const out,
		   const Mat_rm_type type)
{

	int i, j;

	// Resize the output
	out->num_rows = in->num_rows;
	out->num_cols = in->num_cols;
	out->type = type;
	if (resize_Mat_rm(out))
		return SIFT3D_FAILURE;

#define CONVERT_TYPE(type_in, type_out) \
    SIFT3D_MAT_RM_LOOP_START(in, i, j) \
    SIFT3D_MAT_RM_GET(out, i, j, type_out) = (type_out) \
        SIFT3D_MAT_RM_GET(in, i, j, type_in); \
    SIFT3D_MAT_RM_LOOP_END

#define CONVERT_TYPE_OUTER(type_out) \
    switch (in->type) { \
    case SIFT3D_DOUBLE: \
        CONVERT_TYPE(double, type_out) \
        break; \
    case SIFT3D_FLOAT: \
        CONVERT_TYPE(float, type_out) \
        break; \
    case SIFT3D_INT: \
        CONVERT_TYPE(int, type_out) \
        break; \
    default: \
        puts("convert_Mat_rm: unknown type of input matrix \n"); \
        return SIFT3D_FAILURE; \
    }

	// Convert to the specified type
	switch (type) {
	case SIFT3D_DOUBLE:
		CONVERT_TYPE_OUTER(double) break;
	case SIFT3D_FLOAT:
		CONVERT_TYPE_OUTER(float) break;
	case SIFT3D_INT:
		CONVERT_TYPE_OUTER(int) break;
	default:
		puts("convert_Mat_rm: unknown destination type \n");
		return SIFT3D_FAILURE;
	}

#undef CONVERT_TYPE_OUTER
#undef CONVERT_TYPE

	return SIFT3D_SUCCESS;
}

/* Shortcut function to initalize a matrix.
 * 
 * Parameters:
 *      mat - The matrix to be initialized
 *      num_rows - The number of rows
 *      num_cols - The number of columns
 *      type - The data type 
 *      set_zero - If true, initializes the elements to zero.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int init_Mat_rm(Mat_rm *const mat, const int num_rows, const int num_cols,
                const Mat_rm_type type, const int set_zero) {

        mat->type = type;
        mat->num_rows = num_rows;
        mat->num_cols = num_cols;
        mat->u.data_double = NULL;
        mat->size = 0;
        mat->static_mem = SIFT3D_FALSE;

        if (resize_Mat_rm(mat))
                return SIFT3D_FAILURE;
        
        if (set_zero && zero_Mat_rm(mat))
                return SIFT3D_FAILURE;
    
        return SIFT3D_SUCCESS;
}

/* As init_Mat_rm, but aliases data memory with pointer p. The flag 
 * mat->static_mem is set, and the matrix does not need to be freed with 
 * cleanup_Mat_rm. But, an error will be thrown if the user attempts to resize
 * the memory. That is, resize_Mat_rm will only return success if the size of 
 * the matrix does not change. */ 
int init_Mat_rm_p(Mat_rm *const mat, const void *const p, const int num_rows, 
                  const int num_cols, const Mat_rm_type type, 
                  const int set_zero) {

        // Perform normal initialization
        if (init_Mat_rm(mat, num_rows, num_cols, type, set_zero))
                return SIFT3D_FAILURE;

        // Clean up any existing memory
        cleanup_Mat_rm(mat);

        // Alias with provided memory and set the static flag
        mat->u.data_double = (double *) p;
        mat->static_mem = SIFT3D_TRUE;

        // Optionally set to zero 
        if (set_zero && zero_Mat_rm(mat))
                return SIFT3D_FAILURE;

        return SIFT3D_SUCCESS;
}

/* Prints the type of mat into the string str. */
void sprint_type_Mat_rm(const Mat_rm * const mat, char *const str)
{
	switch (mat->type) {
	case SIFT3D_DOUBLE:
		sprintf(str, "double");
		break;
	case SIFT3D_FLOAT:
		sprintf(str, "float");
		break;
	case SIFT3D_INT:
		sprintf(str, "int");
		break;
	default:
		sprintf(str, "<sprint_type_Mat_rm: unknown type>");
	}
}

/* Concatenate two matrices. If dim = 0, concatenates vertically, i.e. 
 *      dst = [src1
 *             src2].
 * If dim = 1, concatenates horizontally, i.e  
 *      dst = [src1 src2]. 
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE on failure. */
int concat_Mat_rm(const Mat_rm * const src1, const Mat_rm * const src2,
		    Mat_rm * const dst, const int dim)
{

        int off[2], dst_dims[2];
	int i, j;

        const int dims1[] = {src1->num_rows, src1->num_cols};
        const int dims2[] = {src2->num_rows, src2->num_cols};

	// Verify inputs
	if (dims1[1 - dim] != dims2[1 - dim]) {
		SIFT3D_ERR("concat_Mat_rm: incompatible dimensions: "
			"left: [%d x %d] right: [%d x %d] dim: %d \n", 
                        dims1[0], dims1[1], dims2[0], dims2[1], dim);
		return SIFT3D_FAILURE;
	}
	if (src1->type != src2->type) {

		char type1[1024], type2[1024];

		sprint_type_Mat_rm(src1, type1);
		sprint_type_Mat_rm(src2, type2);

		SIFT3D_ERR("concat_Mat_rm: incompatible types: "
			"left: <%s> right: <%s> \n", type1, type2);

		return SIFT3D_FAILURE;
	}

        // Compute the destination dimensions
        for (i = 0; i < 2; i++) {
                dst_dims[i] = dim == i ? dims1[i] + dims2[i] : dims1[i];
        }

	// Resize dst
	dst->type = src1->type;
	dst->num_rows = dst_dims[0];
	dst->num_cols = dst_dims[1];
	if (resize_Mat_rm(dst))
		return SIFT3D_FAILURE;

        // Compute the offsets
        for (i = 0; i < 2; i++) {
                off[i] = dim == i ? dims1[i] : 0;
        }

#define COPY_DATA(type) \
        /* Copy src1 data */ \
        SIFT3D_MAT_RM_LOOP_START(src1, i, j) \
                SIFT3D_MAT_RM_GET(dst, i, j, type) = \
                        SIFT3D_MAT_RM_GET(src1, i, j, type); \
        SIFT3D_MAT_RM_LOOP_END \
        \
        /* Copy src2 data */ \
        SIFT3D_MAT_RM_LOOP_START(src2, i, j) \
        \
                SIFT3D_MAT_RM_GET(dst, i + off[0], j + off[1], type) = \
                        SIFT3D_MAT_RM_GET(src2, i, j, type); \
        \
        SIFT3D_MAT_RM_LOOP_END

	// Copy the data
	switch (dst->type) {
	case SIFT3D_DOUBLE:
		COPY_DATA(double);
		break;
	case SIFT3D_FLOAT:
		COPY_DATA(float);
		break;
	case SIFT3D_INT:
		COPY_DATA(int);
		break;
	default:
		SIFT3D_ERR("concat_Mat_rm: unknown type \n");
		return SIFT3D_FAILURE;
	}

#undef COPY_DATA

	return SIFT3D_SUCCESS;
}

/* Copies a matrix. dst will be resized. */
int copy_Mat_rm(const Mat_rm * const src, Mat_rm * const dst)
{

	// Resize dst
	dst->type = src->type;
	dst->num_rows = src->num_rows;
	dst->num_cols = src->num_cols;
	if (resize_Mat_rm(dst))
		return SIFT3D_FAILURE;

	// Copy the data (use memmove because of static mode)
	memmove(dst->u.data_double, src->u.data_double, src->size);

	return SIFT3D_SUCCESS;
}

/* Print a matrix to stdout. The matrix must be initialized. */
int print_Mat_rm(const Mat_rm * const mat)
{

	int i, j;

#define PRINT_MAT_RM(type, format) \
    SIFT3D_MAT_RM_LOOP_START(mat, i, j) \
        printf("%" #format " ", SIFT3D_MAT_RM_GET(mat, i, j, type)); \
        SIFT3D_MAT_RM_LOOP_COL_END \
        puts("\n"); \
    SIFT3D_MAT_RM_LOOP_ROW_END

	switch (mat->type) {
	case SIFT3D_DOUBLE:
		PRINT_MAT_RM(double, f) break;
	case SIFT3D_FLOAT:
		PRINT_MAT_RM(float, f) break;
	case SIFT3D_INT:
		PRINT_MAT_RM(int, d) break;
	default:
		puts("print_Mat_rm: unknown type \n");
		return SIFT3D_FAILURE;
	}
#undef PRINT_MAT_RM

	return SIFT3D_SUCCESS;
}

/* Re-sizes a matrix. The following fields
 * must already be initialized:
 * -num_rows
 * -num_cols
 * -type
 * -u.data_* (NULL for first use, non-null for resize)
 *
 * The following fields will be modified:
 * -size
 * -u.data_* (Change is not guaranteed)
 * 
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise.
 */
int resize_Mat_rm(Mat_rm *const mat) {

    size_t type_size, total_size;

    const int num_rows = mat->num_rows;
    const int num_cols = mat->num_cols;
    double **const data = &mat->u.data_double;
    const size_t numel = num_rows * num_cols;
    const Mat_rm_type type = mat->type;

    // Get the size of the underyling datatype
    switch (type) {
        case SIFT3D_DOUBLE:
            type_size = sizeof(double);
            break;
        case SIFT3D_FLOAT:
            type_size = sizeof(float);
            break;
        case SIFT3D_INT:
            type_size = sizeof(int);
            break;
        default:
            SIFT3D_ERR("resize_Mat_rm: unknown type! \n");
		return SIFT3D_FAILURE;
	}

    // Calculate the new size
    total_size = type_size * numel;

    // Do nothing if the size has not changed
    if (total_size == mat->size)
        return SIFT3D_SUCCESS;
    mat->size = total_size;

    // Check for static reallocation
    if (mat->static_mem) {
        SIFT3D_ERR("resize_Mat_rm: illegal re-allocation of static matrix \n");
        return SIFT3D_FAILURE;
    }

    // Reset if the new size is 0 
    if (total_size == 0) {
	cleanup_Mat_rm(mat);
	return init_Mat_rm(mat, num_rows, num_cols, type, SIFT3D_FALSE);
    }

    // Re-allocate the memory
    if ((*data = (double *) SIFT3D_safe_realloc(*data, total_size)) == NULL) {
        mat->size = 0;
        return SIFT3D_FAILURE;
    }

    return SIFT3D_SUCCESS;
}

/* Set all elements to zero */
int zero_Mat_rm(Mat_rm *const mat)
{

	int i, j;

#define SET_ZERO(type) \
    SIFT3D_MAT_RM_LOOP_START(mat, i, j) \
        SIFT3D_MAT_RM_GET(mat, i, j, type) = (type) 0; \
    SIFT3D_MAT_RM_LOOP_END

	switch (mat->type) {
	case SIFT3D_DOUBLE:
		SET_ZERO(double);
		break;
	case SIFT3D_FLOAT:
		SET_ZERO(float);
		break;
	case SIFT3D_INT:
		SET_ZERO(int);
		break;
	default:
		return SIFT3D_FAILURE;
	}
#undef SET_ZERO

	return SIFT3D_SUCCESS;
}

/* Set a matrix to identity. 
 *
 * Parameters:
 *   n: The length of the square matrix. The output will have size [n x n].
 *   mat: The matrix to be set.
 */
int identity_Mat_rm(const int n, Mat_rm *const mat) {

        int i;

        // Resize the matrix 
        mat->num_rows = mat->num_cols = n;
        if (resize_Mat_rm(mat))
                return SIFT3D_FAILURE;

        // Set to identity
        zero_Mat_rm(mat);

#define SET_IDENTITY(type) \
        for (i = 0; i < n; i++) { \
                SIFT3D_MAT_RM_GET(mat, i, i, type) = (type) 1; \
        }

        SIFT3D_MAT_RM_TYPE_MACRO(mat, identity_Mat_rm_quit, SET_IDENTITY);
#undef SET_IDENTITY

        return SIFT3D_SUCCESS;

identity_Mat_rm_quit:
        return SIFT3D_FAILURE;
}

/* De-allocate the memory for a Mat_rm struct, unless it was initialized in
 * static mode. */
void cleanup_Mat_rm(Mat_rm *mat) {

    if (mat->u.data_double == NULL)
        return;

    if (!mat->static_mem)
        free(mat->u.data_double);
}

/* Make a grid with the specified spacing between lines and line width. 
 * Uses the default stride and initializes grid. */
int draw_grid(Image * grid, int nx, int ny, int nz, int spacing, int line_width)
{

	int x, y, z;

	const double line_half_width = (double)line_width / 2.0;

	// Verify inputs
	if (spacing < 2 || line_width < 1 || line_width > spacing)
		return SIFT3D_FAILURE;

	if (init_im_with_dims(grid, nx, ny, nz, 1))
		return SIFT3D_FAILURE;

	SIFT3D_IM_LOOP_START(grid, x, y, z)
	    if (x % spacing == 0 || y % spacing == 0 || z % spacing == 0) {
		int x_draw, y_draw, z_draw, x_start, x_end, y_start,
		    y_end, z_start, z_end;

		// Draw a line  
		x_start = SIFT3D_MAX(x - line_half_width, 0);
		y_start = SIFT3D_MAX(y - line_half_width, 0);
		z_start = SIFT3D_MAX(z - line_half_width, 0);
		x_end = SIFT3D_MIN(x + line_half_width + 1, nx - 1);
		y_end = SIFT3D_MIN(y + line_half_width + 1, ny - 1);
		z_end = SIFT3D_MIN(z + line_half_width + 1, nz - 1);
		SIFT3D_IM_LOOP_LIMITED_START(grid, x_draw, y_draw, z_draw,
					     x_start, x_end,
					     y_start, y_end, z_start, z_end)
		    if (abs(x_draw - x) < line_half_width &&
			abs(y_draw - y) < line_half_width &&
			abs(z_draw - z) < line_half_width)
			SIFT3D_IM_GET_VOX(grid, x_draw, y_draw,
					  z_draw, 0) = 1.0f;
	SIFT3D_IM_LOOP_END}
	SIFT3D_IM_LOOP_END return SIFT3D_SUCCESS;
}

/* Draw points in in image*/
int draw_points(const Mat_rm * const in, const int *const dims, int radius,
		Image * const out)
{

	Mat_rm in_i;
	int i, x, y, z;

	// Initialize intermediates
	if (init_Mat_rm(&in_i, 0, 0, SIFT3D_INT, SIFT3D_FALSE))
		return SIFT3D_FAILURE;

	// Resize the output
	memcpy(SIFT3D_IM_GET_DIMS(out), dims, IM_NDIMS * sizeof(int));
	out->nc = 1;
	im_default_stride(out);
	if (im_resize(out))
		goto draw_points_quit;
	im_zero(out);

	// Convert the input to integer
	if (convert_Mat_rm(in, &in_i, SIFT3D_INT))
		goto draw_points_quit;

	for (i = 0; i < in->num_rows; i++) {
		const int cx = SIFT3D_MAT_RM_GET(&in_i, i, 0, int);
		const int cy = SIFT3D_MAT_RM_GET(&in_i, i, 1, int);
		const int cz = SIFT3D_MAT_RM_GET(&in_i, i, 2, int);
		const int x_start = SIFT3D_MAX(cx - radius, 0);
		const int y_start = SIFT3D_MAX(cy - radius, 0);
		const int z_start = SIFT3D_MAX(cz - radius, 0);
		const int x_end = SIFT3D_MIN(cx + radius, dims[0] - 1);
		const int y_end = SIFT3D_MIN(cy + radius, dims[1] - 1);
		const int z_end = SIFT3D_MIN(cz + radius, dims[2] - 1);

		// Draw the point
		SIFT3D_IM_LOOP_LIMITED_START(out, x, y, z, x_start, x_end,
					     y_start, y_end, z_start, z_end)
		    SIFT3D_IM_GET_VOX(out, x, y, z, 0) = 1.0f;
	SIFT3D_IM_LOOP_END}

	// Clean up
	cleanup_Mat_rm(&in_i);
	return SIFT3D_SUCCESS;

 draw_points_quit:
	cleanup_Mat_rm(&in_i);
	return SIFT3D_FAILURE;
}

/* Draw lines between two sets of points.
 * TODO currently only does XY plane. Add support for other planes */
int draw_lines(const Mat_rm * const points1, const Mat_rm * const points2,
	       const int *const dims, Image * const out)
{

	Mat_rm points1_d, points2_d;
	double xd;
	int i, y;

	// Parameters
	const double line_step = 0.1;

	// Verify inputs
	if (points1->num_rows != points2->num_rows ||
	    points1->num_cols != points2->num_cols ||
	    points1->num_cols != IM_NDIMS) {
		puts("draw_lines: invalid points dimensions \n");
		return SIFT3D_FAILURE;
	}
	// Initialize intermediates
	if (init_Mat_rm(&points1_d, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE) ||
	    init_Mat_rm(&points2_d, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE))
		return SIFT3D_FAILURE;

	// Resize the output image
	memcpy(SIFT3D_IM_GET_DIMS(out), dims, IM_NDIMS * sizeof(int));
	out->nc = 1;
	im_default_stride(out);
	if (im_resize(out))
		goto draw_lines_quit;
	im_zero(out);

	// Convert the inputs to double
	if (convert_Mat_rm(points1, &points1_d, SIFT3D_DOUBLE) ||
	    convert_Mat_rm(points2, &points2_d, SIFT3D_DOUBLE))
		goto draw_lines_quit;

	for (i = 0; i < points1->num_rows; i++) {

		const double p1x = SIFT3D_MAT_RM_GET(&points1_d, i, 0, double);
		const double p2x = SIFT3D_MAT_RM_GET(&points2_d, i, 0, double);
		const double p1y = SIFT3D_MAT_RM_GET(&points1_d, i, 1, double);
		const double p2y = SIFT3D_MAT_RM_GET(&points2_d, i, 1, double);
		const double p1z = SIFT3D_MAT_RM_GET(&points1_d, i, 2, double);
		const double p2z = SIFT3D_MAT_RM_GET(&points2_d, i, 2, double);

		// Check the bounds
		if (!IM_CONTAINS(out, p1x, p1y, p1z) ||
		    !IM_CONTAINS(out, p2x, p2y, p2z))
			continue;

		// Get the bounds of the line
		const double x_start = SIFT3D_MIN(p1x, p2x) + 0.5;
		const double x_end = SIFT3D_MAX(p1x, p2x) + 0.5;
		const int zi = (int)p1z;

		// Check if the line is vertical
		const int vert = fabs(x_start - x_end) < 1.0;

		// Draw the line
		if (vert) {

			const int xi = (int)x_start;
			const int y_start = (int)SIFT3D_MIN(p1y, p2y);
			const int y_end = (int)SIFT3D_MAX(p1y, p2y);

			for (y = y_start; y <= y_end; y++) {
				SIFT3D_IM_GET_VOX(out, xi, y, zi, 0) = 1.0f;
			}

		} else {

			// Get the line parameters
			const double y_slope = p1x < p2x ? (p2y - p1y) /
			    (p2x - p1x) : (p1y - p2y) / (p1x - p2x);
			const double b = p1y + 0.5 - (p1x + 0.5) * y_slope;

			for (xd = x_start; xd <= x_end; xd += line_step) {

				const double yd = y_slope * xd + b;
				const int xi = (int)xd;
				const int yi = (int)yd;

				if (yi < 0 || yi > dims[1] - 1)
					continue;

				SIFT3D_IM_GET_VOX(out, xi, yi, zi, 0) = 1.0f;
			}
		}
	}

	// Clean up
	cleanup_Mat_rm(&points1_d);
	cleanup_Mat_rm(&points2_d);

	return SIFT3D_SUCCESS;

 draw_lines_quit:
	cleanup_Mat_rm(&points1_d);
	cleanup_Mat_rm(&points2_d);
	return SIFT3D_FAILURE;
}

/* Detect the format of the supplied file name. */
im_format im_get_format(const char *path) {

        struct stat st;
        const char *ext;

        // Check if the file exists and is a directory
        if (stat(path, &st) == 0) {
                if (S_ISDIR(st.st_mode))
                        return DIRECTORY;
        } 

        // If not a directory, get the file extension
        ext = get_file_ext(path);

        // Check the known types
        if (!strcmp(ext, ext_analyze) || !strcmp(ext, ext_gz) ||
                !strcmp(ext, ext_nii))
                return NIFTI;

        if (!strcmp(ext, ext_dcm))
                return DICOM;

        if (!strcmp(ext, ext_dir))
                return DIRECTORY;

        // The type was not recognized
        return UNKNOWN;
}

/* Read an image from a file. The file extension must match one of the
 * supported formats.
 *
 * Supported formats:
 * - Analyze (.img, .img.gz)
 * - DICOM (.dcm)
 * - Directory of DICOM files
 * - NIFTI-1 (.nii, .nii.gz)
 *
 * Return values:
 * -SIFT3D_SUCCESS - Image successfully read
 * -SIFT3D_FILE_DOES_NOT_EXIST - The file does not exist
 * -SIFT3D_UNSUPPORTED_FILE_TYPE - The file type is not supported
 * -SIFT3D_WRAPPER_NOT_COMPILED - The file type is supported, but the wrapper 
 *      library was not compiled.
 * -SIFT3D_UNEVEN_SPACING - The image slices are unevenly spaced.
 * -SIFT3D_INCONSISTENT_AXES  - The image slices have inconsistent axes.
 * -SIFT3D_DUPLICATE_SLICES - There are multiple slices in the same location.
 * -SIFT3D_FAILURE - Other error
 */
int im_read(const char *path, Image *const im) {

        struct stat st;
        int ret;

        // Ensure the file exists
        if (stat(path, &st) != 0) {
                SIFT3D_ERR("im_read: failed to find file %s \n", path);
                return SIFT3D_FILE_DOES_NOT_EXIST;
        }

        // Get the file format and write the file
        switch (im_get_format(path)) {
        case ANALYZE:
        case NIFTI:
                ret = read_nii(path, im);
                break;
        case DICOM:
                ret = read_dcm(path, im);
                break;
        case DIRECTORY:
                ret = read_dcm_dir(path, im);
                break;
        case FILE_ERROR:
                ret = SIFT3D_FAILURE;
                break;
        case UNKNOWN:
        default:
                SIFT3D_ERR("im_read: unrecognized file extension "
                        "from file %s \n", path);
                ret = SIFT3D_UNSUPPORTED_FILE_TYPE;
        }

        return ret;
}

/* Write an image to a file.
 * 
 * Supported formats:
 * -DICOM (.dcm)
 * -Directory of DICOM files
 * -NIFTI (.nii, .nii.gz)
 *
 * Return values:
 * -SIFT3D_SUCCESS - Successfully wrote the image
 * -SIFT3D_UNSUPPORTED_FILE_TYPE - Cannot write this file type
 * -SIFT3D_FAILURE - Other error
 */
int im_write(const char *path, const Image *const im) {

	// Create the path
	if (mkpath(path, out_mode))
		return SIFT3D_FAILURE;

        // Get the file format 
        switch (im_get_format(path)) {
        case ANALYZE:
        case NIFTI:
                return write_nii(path, im);
        case DICOM:
                return write_dcm(path, im, NULL, -1.0f);
        case DIRECTORY:

                // Create the directory
                if (do_mkdir(path, out_mode)) {
                        SIFT3D_ERR("im_write: failed to create directory "
                                "%s \n", path);
                        return SIFT3D_FAILURE;
                }

                return write_dcm_dir(path, im, NULL);
        case UNKNOWN:
        default:
                // Otherwise, the file extension was not found
                SIFT3D_ERR("im_write: unrecognized file extension "
                        "from file %s \n", path);

                return SIFT3D_UNSUPPORTED_FILE_TYPE;
        }

        // Unreachable code
        return SIFT3D_FAILURE;
}

/* Separate the file name component from its path */
static const char *get_file_name(const char *path) {

        const char *name;

        // Get the last file separator
        name = strrchr(path, SIFT3D_FILE_SEP);
        return name == NULL ? path : name;
}

/* Get the extension of a file name */
static const char *get_file_ext(const char *name)
{

        const char *dot;

        // Get the file name component
        name = get_file_name(name);

        // Get the last dot
	dot = strrchr(name, '.');

	return dot == NULL || dot == name ? "" : dot + 1;
}

/* Get the parent directory of a file. The returned string must later be
 * freed. */
char *im_get_parent_dir(const char *path) {

        ptrdiff_t file_pos;
        char *dirName;

        // Duplicate the path so we can edit it
        dirName = strndup(path, FILENAME_MAX); 

        // Subtract away the file name
        file_pos = get_file_name(path) - path;
        if (file_pos > 0)
                dirName[file_pos] = '\0';
       
        return dirName; 
}

/* Write a matrix to a .csv or .csv.gz file. */
int write_Mat_rm(const char *path, const Mat_rm * const mat)
{

	FILE *file;
	gzFile gz;
	const char *ext;
	int i, j, compress;

	const char *mode = "w";

	// Validate and create the output directory
	if (mkpath(path, out_mode))
		return SIFT3D_FAILURE;

	// Get the file extension
	ext = get_file_ext(path);

	// Check if we need to compress the file
	compress = strcmp(ext, ext_gz) == 0;

	// Open the file
	if (compress) {
		if ((gz = gzopen(path, mode)) == Z_NULL)
			return SIFT3D_FAILURE;
	} else {
		if ((file = fopen(path, mode)) == NULL)
			return SIFT3D_FAILURE;
	}

#define WRITE_MAT(mat, format, type) \
    SIFT3D_MAT_RM_LOOP_START(mat, i, j) \
                const char delim = j < mat->num_cols - 1 ? ',' : '\n'; \
                if (compress) { \
                        gzprintf(gz, format, SIFT3D_MAT_RM_GET(mat, i, j, \
                                type)); \
                        gzputc(gz, delim); \
                } else { \
                	fprintf(file, format, SIFT3D_MAT_RM_GET(mat, i, j, \
                 		type)); \
                        fputc(delim, file); \
                } \
    SIFT3D_MAT_RM_LOOP_END

	// Write the matrix
	switch (mat->type) {
	case SIFT3D_DOUBLE:
		WRITE_MAT(mat, "%f", double); 
		break;
	case SIFT3D_FLOAT:
		WRITE_MAT(mat, "%f", float); 
		break;
	case SIFT3D_INT:
		WRITE_MAT(mat, "%d", int); 
		break;
	default:
		goto write_mat_quit;
	}
#undef WRITE_MAT

	// Check for errors and finish writing the matrix
	if (compress) {
		if (gzclose(gz) != Z_OK)
			goto write_mat_quit;
	} else {
		if (ferror(file))
			goto write_mat_quit;
		fclose(file);
	}

	return SIFT3D_SUCCESS;

 write_mat_quit:
	if (compress) {
		gzclose(gz);
	} else {
		fclose(file);
	}
	return SIFT3D_FAILURE;
}

/* Shortcut to initialize an image for first-time use.
 * Allocates memory, and assumes the default stride. This
 * function calls init_im and initializes all values to 0. */
int init_im_with_dims(Image *const im, const int nx, const int ny, const int nz,
        const int nc)
{

	init_im(im);
	im->nx = nx;
	im->ny = ny;
	im->nz = nz;
	im->nc = nc;

	im_default_stride(im);
	if (im_resize(im))
		return SIFT3D_FAILURE;

	im_zero(im);

	return SIFT3D_SUCCESS;
}

/* Calculate the strides of an image object in the default
 * manner. The following parameters must be initialized:
 * -nx
 * -ny
 * -nz
 * -nc
 * If a dimension is not used, its size should be set
 * to 1. */
void im_default_stride(Image *const im)
{

        size_t prod;
	int i;

	prod = (size_t) im->nc;
	SIFT3D_IM_GET_STRIDES(im)[0] = prod;

	for (i = 1; i < IM_NDIMS; i++) {
		prod *= SIFT3D_IM_GET_DIMS(im)[i - 1];
		SIFT3D_IM_GET_STRIDES(im)[i] = prod;
	}
}

/* Pads an image to a new size. Prior to calling this function, initialize 
 * pad with all dimensions and strides, as in im_resize. Other metadata,
 * such as units, will be copied from im to pad. */
int im_pad(const Image * const im, Image * const pad)
{

	int x, y, z, c;

	const int pad_x_end = pad->nx - 1;
	const int pad_y_end = pad->ny - 1;
	const int pad_z_end = pad->nz - 1;

	const int data_x_end = SIFT3D_MIN(im->nx - 1, pad_x_end);
	const int data_y_end = SIFT3D_MIN(im->ny - 1, pad_y_end);
	const int data_z_end = SIFT3D_MIN(im->nz - 1, pad_z_end);

        // Copy relevant metadata, omitting dimensions and strides
        memcpy(SIFT3D_IM_GET_UNITS(pad), SIFT3D_IM_GET_UNITS(im), 
                IM_NDIMS * sizeof(double));

	// Resize the output 
	if (im_resize(pad))
		return SIFT3D_FAILURE;

	// Copy the image 
	SIFT3D_IM_LOOP_LIMITED_START_C(im, x, y, z, c, 0, data_x_end, 0,
				       data_y_end, 0, data_z_end)

	    SIFT3D_IM_GET_VOX(pad, x, y, z, c) =
	    SIFT3D_IM_GET_VOX(im, x, y, z, c);

	SIFT3D_IM_LOOP_END_C
	    // Pad the remaining data with zeros
	    SIFT3D_IM_LOOP_LIMITED_START_C(im, x, y, z, c, data_x_end,
					   pad_x_end, data_y_end, pad_y_end,
					   data_z_end, pad_z_end)

	    SIFT3D_IM_GET_VOX(pad, x, y, z, c) = 0.0f;

	SIFT3D_IM_LOOP_END_C return SIFT3D_SUCCESS;
}

/* Resize an image according to the current nx, ny,
 * and nz. Does not modify scale space information or
 * strides. Prior to calling this function, use init_im(im)
 * and initialize the following fields:
 * -nx
 * -ny
 * -nz
 * -nc
 * -xs (can be set by im_default_stride(im)) 
 * -ys (can be set by im_default_stride(im)) 
 * -zs (can be set by im_default_stride(im)) 
 *
 * All of this initialization can also be done with
 * init_im_with_dims(), which calls this function.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise.
 */
int im_resize(Image *const im)
{

	int i;

	//FIXME: This will not work for strange strides
	const size_t size = im->nx * im->ny * im->nz * im->nc;

	// Verify inputs
	for (i = 0; i < IM_NDIMS; i++) {

		const int dim = SIFT3D_IM_GET_DIMS(im)[i];

		if (dim > 0)
			continue;

		SIFT3D_ERR("im_resize: invalid dimension %d: %d \n", i,
			dim);
		return SIFT3D_FAILURE;
	}
	if (im->nc < 1) {
		SIFT3D_ERR("im_resize: invalid number of channels: %d \n",
			im->nc);
		return SIFT3D_FAILURE;
	}

        // Do nothing if the size has not changed
        if (im->size == size)
                return SIFT3D_SUCCESS;
	im->size = size;

	// Allocate new memory
	im->data = SIFT3D_safe_realloc(im->data, size * sizeof(float));

#ifdef SIFT3D_USE_OPENCL
	{
		cl_int err;
		int initialized;

		if (cl_data.valid) {
			initialized = (im->data != NULL);

			// Destroy the old image
			if (initialized && im->cl_valid)
				clReleaseMemObject(im->cl_image);

			// Init an OpenCL image
			if (im->nz > 0) {
				im->cl_image = clCreateImage2D(cl_data.context,
							       cl_data.
							       mem_flags,
							       &cl_data.
							       image_format,
							       im->nx, im->ny,
							       im->ys,
							       im->data, &err);
			} else {
				im->cl_image = clCreateImage3D(cl_data.context,
							       cl_data.
							       mem_flags,
							       &cl_data.
							       image_format,
							       im->nx, im->ny,
							       im->nz,
							       im->ys,
							       im->zs,
							       im->data, &err);
			}

			if (err != CL_SUCCESS) {
				im->cl_valid = SIFT3D_FALSE;
				return SIFT3D_FAILURE;
			}

			im->cl_valid = SIFT3D_TRUE;
		}
	}
#endif
	return size != 0 && im->data == NULL ? SIFT3D_FAILURE : SIFT3D_SUCCESS;
}

/* Concatenate two images in dimension dim, so that src1 comes before src2.
 * For example, if dim == 0, the images are horizontally concatenated in x,
 * so that src1 is on the left and src2 is on the right. Resizes dst. 
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int im_concat(const Image * const src1, const Image * const src2, const int dim,
	      Image * const dst)
{

	int off[IM_NDIMS], dims_out[IM_NDIMS];
	int i, x, y, z, c;

	const int nc = src1->nc;

	// Verify inputs
	for (i = 0; i < IM_NDIMS; i++) {

		const int src1d = SIFT3D_IM_GET_DIMS(src1)[i];
		const int src2d = SIFT3D_IM_GET_DIMS(src2)[i];

		if (i == dim)
			continue;
		if (src1d != src2d) {
			SIFT3D_ERR("im_concat: dimension %d must be "
				"equal in input images. src1: %d src2: %d \n",
				i, src1d, src2d);
			return SIFT3D_FAILURE;
		}
	}
	if (src1->nc != src2->nc) {
		SIFT3D_ERR("im_concat: images must have an equal number "
			"of channels. src1: %d src2: %d \n", src1->nc,
			src2->nc);
		return SIFT3D_FAILURE;
	}
	// Get the output dimensions and offsets
	for (i = 0; i < IM_NDIMS; i++) {

		const int src1d = SIFT3D_IM_GET_DIMS(src1)[i];
		const int src2d = SIFT3D_IM_GET_DIMS(src2)[i];

                dims_out[i] = dim == i ? src1d + src2d : src1d;
                off[i] = dim == i ? src1d : 0;
	}

	// Resize dst 
	memcpy(SIFT3D_IM_GET_DIMS(dst), dims_out, IM_NDIMS * sizeof(int));
	dst->nc = nc;
	im_default_stride(dst);
	if (im_resize(dst))
		return SIFT3D_FAILURE;

	// Copy the data from src1 
	SIFT3D_IM_LOOP_START_C(src1, x, y, z, c)
	    SIFT3D_IM_GET_VOX(dst, x, y, z, c) =
	    SIFT3D_IM_GET_VOX(src1, x, y, z, c);
	SIFT3D_IM_LOOP_END_C 

        // Copy the data from src2
        SIFT3D_IM_LOOP_START_C(src2, x, y, z, c)

                // Get the destination coordinates with offsets
                const int x_dst = x + off[0];
                const int y_dst = y + off[1];
                const int z_dst = z + off[2];

                // Copy the data from src2
                SIFT3D_IM_GET_VOX(dst, x_dst, y_dst, z_dst, c) =
                    SIFT3D_IM_GET_VOX(src2, x, y, z, c);

        SIFT3D_IM_LOOP_END_C 
        
        return SIFT3D_SUCCESS;
}

/* Upsample an image by a factor of 2 in each dimension.
 * This function resizes dst and modifies its units. */
int im_upsample_2x(const Image *const src, Image *const dst)
{

	double units[IM_NDIMS];
	int dims[IM_NDIMS];
	int i, x, y, z, c, sx, sy, sz;

	const int nc = src->nc;
	const int w = 2;
	const float weight = (float) (1.0 / pow((double) w, IM_NDIMS));

	// Resize the output 
	for (i = 0; i < IM_NDIMS; i++) {
		dims[i] = SIFT3D_IM_GET_DIMS(src)[i] * 2;	
		units[i] = SIFT3D_IM_GET_UNITS(src)[i] / 2.0;
	}
	memcpy(SIFT3D_IM_GET_DIMS(dst), dims, IM_NDIMS * sizeof(int));
	memcpy(SIFT3D_IM_GET_UNITS(dst), units, IM_NDIMS * sizeof(float));
	dst->nc = nc;
	im_default_stride(dst);
	if (im_resize(dst))
		return SIFT3D_FAILURE;

	// TODO: 3-pass (separable) upsample might be faster

	// Upsample
	SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)

		const int sx_start = x >> 1;
		const int sy_start = y >> 1;
		const int sz_start = z >> 1;

		const int sx_end = sx_start + w - 1;
		const int sy_end = sy_start + w - 1;
		const int sz_end = sz_start + w - 1;

		SIFT3D_IM_GET_VOX(dst, x, y, z, c) = 0;

		SIFT3D_IM_LOOP_LIMITED_START(dst, sx, sy, sz, sx_start,
					       sx_end, sy_start, sy_end,
					       sz_start, sz_end)

			SIFT3D_IM_GET_VOX(dst, x, y, z, c) += 
				SIFT3D_IM_GET_VOX(src, sx, sy, sz, c);

		SIFT3D_IM_LOOP_END 

		SIFT3D_IM_GET_VOX(dst, x, y, z, c) *= weight;	

	SIFT3D_IM_LOOP_END_C

	return SIFT3D_SUCCESS;
}

/* Downsample an image by a factor of 2 in each dimension.
 * This function initializes dst with the proper 
 * dimensions, and allocates memory. */
int im_downsample_2x(const Image *const src, Image *const dst)
{

	int x, y, z, c;

	// Initialize dst
	dst->nx = (int)floor((double)src->nx / 2.0);
	dst->ny = (int)floor((double)src->ny / 2.0);
	dst->nz = (int)floor((double)src->nz / 2.0);
	dst->nc = src->nc;
	im_default_stride(dst);
	if (im_resize(dst))
		return SIFT3D_FAILURE;

	// Downsample
	SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)

		const int src_x = x << 1;
		const int src_y = y << 1;
		const int src_z = z << 1;

		SIFT3D_IM_GET_VOX(dst, x, y, z, c) =
		    SIFT3D_IM_GET_VOX(src, src_x, src_y, src_z, c);
	SIFT3D_IM_LOOP_END_C 

	return SIFT3D_SUCCESS;
}

/* Same as im_downsample_2x, but with OpenCL acceleration. This function DOES NOT
 * read the results back dst->data. Use im_read_back for that. */
int im_downsample_2x_cl(Image * src, Image * dst)
{
#ifdef SIFT3D_USE_OPENCL
	size_t global_work_size[3];
	cl_int err, dim;
	cl_kernel kernel;

	// Verify image dimensions
	if (src->nx % 2 != 0 || src->ny % 2 != 0 || src->nz % 2 != 0)
		return SIFT3D_FAILURE;

	// Initialize dst dimensions, resized in im_set_kernel_arg
	dst->nx = src->nx / 2;
	dst->ny = src->ny / 2;
	dst->nz = src->nz / 2;
	dst->nc = src->nc;
	im_default_stride(dst);

	// Do not have a 2D kernel right now
	assert(src->nz > 0);
	dim = 3;
	global_work_size[0] = dst->nx;
	global_work_size[1] = dst->ny;
	global_work_size[2] = dst->nz;

	kernel = cl_data.kernels.downsample_2x_3d;
	im_set_kernel_arg(kernel, 0, src);
	im_set_kernel_arg(kernel, 1, dst);

	err =
	    clEnqueueNDRangeKernel(cl_data.queues[0], kernel, dim, NULL,
				   global_work_size, NULL, 0, NULL, NULL);
	return (int)err;
#else
	printf
	    ("im_downsample_2x_cl: This version was not compiled with OpenCL!\n");
	return SIFT3D_FAILURE;
#endif
}

/* Loads the C-accessible data of an image into its OpenCL data. If blocking is set,
 * block the function until the load is complete. */
int im_load_cl(Image * im, int blocking)
{
#ifdef SIFT3D_USE_OPENCL
	const size_t origin[] = { 0, 0, 0 };
	const size_t region[] = { im->nx, im->ny, im->nz };
	const cl_bool cl_blocking = (blocking) ? CL_TRUE : CL_FALSE;
	return clEnqueueWriteImage(cl_data.queues[0], im->cl_image, cl_blocking,
				   origin, region, im->ys, im->zs,
				   im->data, 0, NULL, NULL);
#else
	printf("im_load_cl: This version was not compiled with OpenCL!\n");
	return SIFT3D_FAILURE;
#endif
}

/* Reads the OpenCL data of an image back to its C-accessible data. If blocking is set,
 * block the function until the read is complete. */
int im_read_back(Image * im, int blocking)
{
#ifdef SIFT3D_USE_OPENCL
	const size_t origin[] = { 0, 0, 0 };
	const size_t region[] = { im->nx, im->ny, im->nz };
	const cl_bool cl_blocking = (blocking) ? CL_TRUE : CL_FALSE;
	return clEnqueueReadImage(cl_data.queues[0], im->cl_image, cl_blocking,
				  origin, region, im->ys, im->zs,
				  im->data, 0, NULL, NULL);
#else
	printf("im_read_back: This version was not compiled with OpenCL!\n");
	return SIFT3D_FAILURE;
#endif
}

/* Updates an Image struct's OpenCL data, if neccessary, and sets it as argument n
 * in the provided kernel. */
int im_set_kernel_arg(cl_kernel kernel, int n, Image * im)
{
#ifdef SIFT3D_USE_OPENCL
	cl_int err;

	if (!im->cl_valid && im_resize(im))
		return SIFT3D_FAILURE;
	err = clSetKernelArg(kernel, n, sizeof(cl_mem), &im->cl_image);
	check_cl_error(err, "im_set_kernel_arg");

	return SIFT3D_SUCCESS;
#else
	printf
	    ("im_set_kernel_arg: This version was not compiled with OpenCL!\n");
	return SIFT3D_FAILURE;
#endif
}

/* Copy an image's dimensions and stride into another. 
 * This function resizes dst.
 * 
 * @param src The source image.
 * @param dst The destination image.
 * @return Returns SIFT3D_SUCCESS or SIFT3D_FAILURE.
 */
int im_copy_dims(const Image * const src, Image * dst)
{
        if (src->data == NULL)
                return SIFT3D_FAILURE;

	dst->nx = src->nx;
	dst->ny = src->ny;
	dst->nz = src->nz;
	dst->xs = src->xs;
	dst->ys = src->ys;
	dst->zs = src->zs;
	dst->nc = src->nc;
        dst->ux = src->ux;
        dst->uy = src->uy;
        dst->uz = src->uz;

	return im_resize(dst);
}

/* Copy an image's data into another. This function
 * changes the dimensions and stride of dst,
 * and allocates memory. */
int im_copy_data(const Image * const src, Image * const dst)
{

	int x, y, z, c;

        // Return if src has no data 
        if (src->data == NULL)
                return SIFT3D_FAILURE;

	// Return if src and dst are the same
	if (dst->data == src->data)
		return SIFT3D_SUCCESS;

	// Resize dst
	if (im_copy_dims(src, dst))
		return SIFT3D_FAILURE;

	// Copy data
        SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)
                SIFT3D_IM_GET_VOX(dst, x, y, z, c) = 
                        SIFT3D_IM_GET_VOX(src, x, y, z, c);
        SIFT3D_IM_LOOP_END_C

	return SIFT3D_SUCCESS;
}

/* Clean up memory for an Image */
void im_free(Image * im)
{
	if (im->data != NULL)
		free(im->data);
}

/* Make a deep copy of a single channel of an image. */
int im_channel(const Image * const src, Image * const dst,
	       const unsigned int chan)
{

	int x, y, z;

	const int c = chan;

	// Verify inputs
	if (c >= src->nc) {
		SIFT3D_ERR("im_channel: invalid channel: %d, image has "
			"%d channels", c, src->nc);
		return SIFT3D_FAILURE;
	}
	// Resize the output
	memcpy(SIFT3D_IM_GET_DIMS(dst), SIFT3D_IM_GET_DIMS(src), 
                IM_NDIMS * sizeof(int));
	dst->nc = 1;
	im_default_stride(dst);
	if (im_resize(dst))
		return SIFT3D_FAILURE;

	// Copy the channel
	SIFT3D_IM_LOOP_START(dst, x, y, z)
	    SIFT3D_IM_GET_VOX(dst, x, y, z, 0) =
	    SIFT3D_IM_GET_VOX(src, x, y, z, c);
	SIFT3D_IM_LOOP_END return SIFT3D_SUCCESS;
}

/* Find the maximum absolute value of an image */
float im_max_abs(const Image *const im) {

        float max;
        int x, y, z, c;

	max = 0.0f;
	SIFT3D_IM_LOOP_START_C(im, x, y, z, c)

	        const float samp = fabsf(SIFT3D_IM_GET_VOX(im, x, y, z, c));
                max = SIFT3D_MAX(max, samp);

	SIFT3D_IM_LOOP_END_C

        return max;
}

/* Scale an image to the [-1, 1] range, where
 * the largest absolute value is 1. */
void im_scale(const Image *const im)
{

	int x, y, z, c;

        // Find the maximum absolute value
	const float max = im_max_abs(im);
        if (max == 0.0f)
	        return;

	// Divide by the max 
	SIFT3D_IM_LOOP_START_C(im, x, y, z, c)
	        SIFT3D_IM_GET_VOX(im, x, y, z, c) /= max;
        SIFT3D_IM_LOOP_END_C
}

/* Subtract src2 from src1, saving the result in
 * dst.
 * Resizes dst. 
 */
int im_subtract(Image * src1, Image * src2, Image * dst)
{

	int x, y, z, c;

	// Verify inputs
	if (src1->nx != src2->nx ||
	    src1->ny != src2->ny ||
	    src1->nz != src2->nz || src1->nc != src2->nc)
		return SIFT3D_FAILURE;

	// Resize the output image
	if (im_copy_dims(src1, dst))
		return SIFT3D_FAILURE;

	SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)
	    SIFT3D_IM_GET_VOX(dst, x, y, z, c) =
	    SIFT3D_IM_GET_VOX(src1, x, y, z, c) -
	    SIFT3D_IM_GET_VOX(src2, x, y, z, c);
	SIFT3D_IM_LOOP_END_C return SIFT3D_SUCCESS;
}

/* Zero an image. */
void im_zero(Image * im)
{

	int x, y, z, c;

	SIFT3D_IM_LOOP_START_C(im, x, y, z, c)
	    SIFT3D_IM_GET_VOX(im, x, y, z, c) = 0.0f;
SIFT3D_IM_LOOP_END_C}

/* Transform an image according to the inverse of the provided tform. 
 * 
 * Paramters:
 *   tform: The transformation. 
 *   src: The input image.
 *   interp: The type of interpolation.
 *   resize: If true, resizes the dst to be the same size as src. Otherwise,
 *     uses the dimensions of dst. 
 *   dst: The output image.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int im_inv_transform(const void *const tform, const Image * const src,
		     const interp_type interp, const int resize, 
                     Image *const dst)
{
	int x, y, z, c;

	// Optionally resize the output image
	if (resize && im_copy_dims(src, dst))
		return SIFT3D_FAILURE;

#define IMUTIL_RESAMPLE(arg) \
    SIFT3D_IM_LOOP_START(dst, x, y, z) \
\
	double transx, transy, transz; \
\
        apply_tform_xyz(tform, (double)x, (double)y, (double)z, \
            &transx, &transy, &transz); \
                        \
                for (c = 0; c < dst->nc; c++) { \
                SIFT3D_IM_GET_VOX(dst, x, y, z, c) = resample_ ## arg(src, \
                                transx, transy, transz, c); \
                } \
    SIFT3D_IM_LOOP_END

	// Transform
	switch (interp) {
	case LINEAR:
		IMUTIL_RESAMPLE(linear)
		    break;
	case LANCZOS2:
		IMUTIL_RESAMPLE(lanczos2)
		    break;
	default:
		SIFT3D_ERR("im_inv_transform: unrecognized "
			"interpolation type");
		return SIFT3D_FAILURE;
	}

#undef RESAMPLE

	return SIFT3D_SUCCESS;
}

/* Helper routine for image transformation. Performs trilinear
 * interpolation, setting out-of-bounds voxels to zero. */
static double resample_linear(const Image * const in, const double x,
			      const double y, const double z, const int c)
{

	// Detect out-of-bounds
	if (x < 0 || x > in->nx - 1 ||
	    y < 0 || y > in->ny - 1 || z < 0 || z > in->nz - 1)
		return 0.0;

	int fx = (int)floor(x);
	int fy = (int)floor(y);
	int fz = (int)floor(z);
	int cx = (int)ceil(x);
	int cy = (int)ceil(y);
	int cz = (int)ceil(z);

	double dist_x = x - fx;
	double dist_y = y - fy;
	double dist_z = z - fz;

	double c0 = SIFT3D_IM_GET_VOX(in, fx, fy, fz, c);
	double c1 = SIFT3D_IM_GET_VOX(in, fx, cy, fz, c);
	double c2 = SIFT3D_IM_GET_VOX(in, cx, fy, fz, c);
	double c3 = SIFT3D_IM_GET_VOX(in, cx, cy, fz, c);
	double c4 = SIFT3D_IM_GET_VOX(in, fx, fy, cz, c);
	double c5 = SIFT3D_IM_GET_VOX(in, fx, cy, cz, c);
	double c6 = SIFT3D_IM_GET_VOX(in, cx, fy, cz, c);
	double c7 = SIFT3D_IM_GET_VOX(in, cx, cy, cz, c);

	double out = c0 * (1.0 - dist_x) * (1.0 - dist_y) * (1.0 - dist_z)
	    + c1 * (1.0 - dist_x) * dist_y * (1.0 - dist_z)
	    + c2 * dist_x * (1.0 - dist_y) * (1.0 - dist_z)
	    + c3 * dist_x * dist_y * (1.0 - dist_z)
	    + c4 * (1.0 - dist_x) * (1.0 - dist_y) * dist_z
	    + c5 * (1.0 - dist_x) * dist_y * dist_z
	    + c6 * dist_x * (1.0 - dist_y) * dist_z
	    + c7 * dist_x * dist_y * dist_z;

	return out;
}

/* Helper routine to resample an image at a point, using the Lanczos kernel */
static double resample_lanczos2(const Image * const im, const double x,
				const double y, const double z, const int c)
{

	double val;
	int xs, ys, zs;

	//TODO: faster separable implementation

	// Kernel parameter
	const double a = 2;

	// Check bounds
	const double xMin = 0;
	const double yMin = 0;
	const double zMin = 0;
	const double xMax = im->nx - 1;
	const double yMax = im->ny - 1;
	const double zMax = im->nz - 1;
	if (x < xMin || y < yMin || z < zMin ||
	    x > xMax || y > yMax || z > zMax)
		return 0.0;

	// Window 
	const int x_start = SIFT3D_MAX(floor(x) - a, xMin);
	const int x_end = SIFT3D_MIN(floor(x) + a, xMax);
	const int y_start = SIFT3D_MAX(floor(y) - a, yMin);
	const int y_end = SIFT3D_MIN(floor(y) + a, yMax);
	const int z_start = SIFT3D_MAX(floor(z) - a, zMin);
	const int z_end = SIFT3D_MIN(floor(z) + a, zMax);

	// Iterate through the window 
	val = 0.0;
	SIFT3D_IM_LOOP_LIMITED_START(in, xs, ys, zs, x_start, x_end, y_start,
				     y_end, z_start, z_end)

        // Evalutate the kernel
	const double xw = fabs((double)xs - x) + DBL_EPSILON;
	const double yw = fabs((double)ys - y) + DBL_EPSILON;
	const double zw = fabs((double)zs - z) + DBL_EPSILON;
	const double kernel = lanczos(xw, a) * lanczos(yw, a) * lanczos(zw, a);

	// Accumulate
	val += kernel * SIFT3D_IM_GET_VOX(im, xs, ys, zs, c);

	SIFT3D_IM_LOOP_END return val;
}

/* Lanczos kernel function */
static double lanczos(double x, double a)
{
	const double pi_x = M_PI * x;
	return a * sin(pi_x) * sin(pi_x / a) / (pi_x * pi_x);
}

/* Resample an image to different units.
 *
 * Parameters:
 *   src: The input image.
 *   units: The new units. 
 *   interp: The type of interpolation to use.
 *   dst: The output image. 
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int im_resample(const Image *const src, const double *const units, 
	const interp_type interp, Image *const dst) {

	Affine aff;
	Mat_rm A;
	double factors[IM_NDIMS];
	int i;

	// Initialize intermediates
	if (init_Mat_rm(&A, IM_NDIMS, IM_NDIMS + 1, SIFT3D_DOUBLE, 
		SIFT3D_TRUE) ||
		init_Affine(&aff, IM_NDIMS))
		return SIFT3D_FAILURE;

	// Compute the scaling factors
	for (i = 0; i < IM_NDIMS; i++) {
		factors[i] = SIFT3D_IM_GET_UNITS(src)[i] / units[i];
	}

	// Set the transformation matrix
	for (i = 0; i < IM_NDIMS; i++) {
		SIFT3D_MAT_RM_GET(&A, i, i, double) = 1.0 / factors[i];
	}
	if (Affine_set_mat(&A, &aff))
		goto im_resample_quit;

	// Set the output dimensions
	dst->nc = src->nc;
	for (i = 0; i < IM_NDIMS; i++) {
		SIFT3D_IM_GET_DIMS(dst)[i] = (int) ceil((double) 
                        SIFT3D_IM_GET_DIMS(src)[i] * factors[i]);
	}
	im_default_stride(dst);	
        if (im_resize(dst))
		goto im_resample_quit;

	// Apply the transformation
	if (im_inv_transform(&aff, src, interp, SIFT3D_FALSE, dst))
		goto im_resample_quit;

        // Set the new output units
	memcpy(SIFT3D_IM_GET_UNITS(dst), units, IM_NDIMS * sizeof(double));

	// Clean up
	cleanup_tform(&aff);
	cleanup_Mat_rm(&A);

	return SIFT3D_SUCCESS;

im_resample_quit:
	cleanup_tform(&aff);
	cleanup_Mat_rm(&A);
	return SIFT3D_FAILURE;
}

/* Horizontally convolves a separable filter with an image, 
 * on CPU. Currently only works in 3D.
 * 
 * This function chooses among the best variant of convolve_sep* based on
 * compilation options and filter parameters.
 * 
 * Parameters: 
 * src - input image (initialized)
 * dst - output image (initialized) 
    int x, y, z;
 * f - filter to be applied
 * dim - dimension in which to convolve
 * unit - the spacing of the filter coefficients
 */
static int convolve_sep(const Image * const src,
			Image * const dst, const Sep_FIR_filter * const f,
			const int dim, const double unit) {

#ifdef SIFT3D_USE_OPENCL
        return convolve_sep_cl(src, dst, f, dim, unit);
#else
	return f->symmetric ? 
                convolve_sep_sym(src, dst, f, dim, unit) : 
                convolve_sep_gen(src, dst, f, dim, unit);
#endif
}

/* Convolve_sep for general filters */
static int convolve_sep_gen(const Image * const src,
			Image * const dst, const Sep_FIR_filter * const f,
			const int dim, const double unit)
{
	register int x, y, z, c, d;

	register const int half_width = f->width / 2;
	register const int nx = src->nx;
	register const int ny = src->ny;
	register const int nz = src->nz;
        register const float conv_eps = 0.1f;
	register const int dim_end = SIFT3D_IM_GET_DIMS(src)[dim] - 1;
        register const float unit_factor =  unit /
                SIFT3D_IM_GET_UNITS(src)[dim];
        register const int unit_half_width = 
                (int) ceilf(half_width * unit_factor);
        int start[] = {0, 0, 0};
        int end[] = {nx - 1, ny - 1, nz - 1};

        // Compute starting and ending points for the convolution dimension
        start[dim] += unit_half_width;
        end[dim] -= unit_half_width + 1;

	//TODO: Convert this to convolve_x, which only convolves in x,
	// then make a wrapper to restride, transpose, convolve x, and transpose 
	// back

	// Resize the output, with the default stride
        if (im_copy_dims(src, dst))
                return SIFT3D_FAILURE;
	im_default_stride(dst);
	if (im_resize(dst))
		return SIFT3D_FAILURE;

	// Initialize the output to zeros
	im_zero(dst);

#define SAMP_AND_ACC(src, dst, tap, coords, c) \
{ \
        float frac; \
\
        const int idx_lo[] = {(coords)[0], (coords)[1], (coords)[2]}; \
        int idx_hi[] = {idx_lo[0], idx_lo[1], idx_lo[2]}; \
\
        /* Convert the physical coordinates to integer indices*/ \
        idx_hi[dim] += 1; \
        frac = (coords)[dim] - (float) idx_lo[dim]; \
\
        /* Sample with linear interpolation */ \
        SIFT3D_IM_GET_VOX(dst, x, y, z, c) += (tap) * \
                ((1.0f - frac) * \
                SIFT3D_IM_GET_VOX(src, idx_lo[0], idx_lo[1], idx_lo[2], c) + \
                frac * \
                SIFT3D_IM_GET_VOX(src, idx_hi[0], idx_hi[1], idx_hi[2], c)); \
}

	// First pass: process the interior
#pragma omp parallel for private(x) private(y) private(c)
	SIFT3D_IM_LOOP_LIMITED_START_C(dst, x, y, z, c, start[0], end[0], 
                start[1], end[1], start[2], end[2])

                float coords[] = { x, y, z };

                for (d = -half_width; d <= half_width; d++) {

                        const float tap = f->kernel[d + half_width];
                        const float step = d * unit_factor;

                        // Adjust the sampling coordinates
                        coords[dim] -= step;

                        // Sample
                        SAMP_AND_ACC(src, dst, tap, coords, c);

                        // Reset the sampling coordinates
                        coords[dim] += step;
                }

	SIFT3D_IM_LOOP_END_C

        // Second pass: process the boundaries
#pragma omp parallel for private(x) private(y) private(c)
        SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)

                const int i_coords[] = { x, y, z };

                // Skip pixels we have already processed
                if (i_coords[dim] >= start[dim] && i_coords[dim] <= end[dim]) 
                        continue;

                // Process the boundary pixel
                for (d = -half_width; d <= half_width; d++) {

                        float coords[] = { x, y, z };
                        const float tap = f->kernel[d + half_width];
                        const float step = d * unit_factor;

                        // Adjust the sampling coordinates
                        coords[dim] -= step;

                        // Mirror coordinates
                        if ((int) coords[dim] < 0) {
                                coords[dim] = -coords[dim];
                                assert((int) coords[dim] >= 0);
                        } else if ((int) coords[dim] >= dim_end) {
                                coords[dim] = 2.0f * dim_end - coords[dim] -    
                                        conv_eps;
                                assert((int) coords[dim] < dim_end);
                        }

                        // Sample
                        SAMP_AND_ACC(src, dst, tap, coords, c);
                }

	SIFT3D_IM_LOOP_END_C 

#undef SAMP_AND_ACC

        return SIFT3D_SUCCESS;
}

/* Same as convolve_sep, but with OpenCL acceleration. This does NOT
 * read back the results to C-accessible data. Use im_read_back for that. */
SIFT3D_IGNORE_UNUSED
static int convolve_sep_cl(const Image * const src, Image * const dst,
			   const Sep_FIR_filter * const f, const int dim,
                           const double unit)
{
#ifdef SIFT3D_USE_OPENCL
	cl_kernel kernel;
	cl_int dx, dy, dz, err;

	const size_t global_work_size[] = { src->nx, src->ny, src->nz };

	// Resize the output, with the default stride
        if (im_copy_dims(src, dst))
                return SIFT3D_FAILURE;
	im_default_stride(dst);
	if (im_resize(dst))
		return SIFT3D_FAILURE;

	// Do not have a 2D kernel right now
	if (dim != 3) {
		printf("convolve_sep_cl: unsupported dimension: %d \n", dim);
	return SIFT3D_FAILURE}

	// Form the dimension offsets
	dx = dy = dz = 0;
	switch (dim) {
	case 0:
		dx = 1;
		break;
	case 1:
		dy = 1;
		break;
	case 2:
		dz = 1;
		break;
	default:
		return SIFT3D_FAILURE;
	}

	kernel = f->cl_apply_unrolled;
	im_set_kernel_arg(kernel, 0, src);
	im_set_kernel_arg(kernel, 1, dst);
	err = clSetKernelArg(kernel, 2, sizeof(cl_int), &dx);
	err |= clSetKernelArg(kernel, 3, sizeof(cl_int), &dy);
	err |= clSetKernelArg(kernel, 4, sizeof(cl_int), &dz);
	check_cl_error(err, "convolve_sep_cl: set kernel arg");

	err =
	    clEnqueueNDRangeKernel(cl_data.queues[0], kernel, dim, NULL,
				   global_work_size, NULL, 0, NULL, NULL);
	return (int)err;
#else
	printf("colvolve_sep_cl: This version was not compiled with OpenCL!\n");
	return SIFT3D_FAILURE;
#endif
}

 /* Convolve_sep for symmetric filters. */
static int convolve_sep_sym(const Image * const src, Image * const dst,
			    const Sep_FIR_filter * const f, const int dim,
                            const double unit)
{

	// TODO: Symmetry-specific function
	return convolve_sep_gen(src, dst, f, dim, unit);
}

/* Permute the dimensions of an image.
 *
 * Arguments: 
 * src - input image (initialized)
 * dim1 - input permutation dimension (x = 0, y = 1, z = 2)
 * dim2 - output permutation dimension (x = 0, y = 1, z = 2)
 * dst - output image (initialized)
 * 
 * example:
 * im_permute(src, dst, 0, 1) -- permute x with y in src
 *                              and save to dst
 */
int im_permute(const Image * const src, const int dim1, const int dim2,
		 Image * const dst)
{
	register int x, y, z, c;

	// Verify inputs
	if (dim1 < 0 || dim2 < 0 || dim1 > 3 || dim2 > 3) {
		printf("im_permute: invalid dimensions: dim1 %d dim2 %d \n",
		       dim1, dim2);
		return SIFT3D_FAILURE;
	}

	// Check for the trivial case
	if (dim1 == dim2) {
		return im_copy_data(src, dst);
        }

        // Permute the units
        memcpy(SIFT3D_IM_GET_UNITS(dst), SIFT3D_IM_GET_UNITS(src), 
                IM_NDIMS * sizeof(double));
        SIFT3D_IM_GET_UNITS(dst)[dim1] = SIFT3D_IM_GET_UNITS(src)[dim2];
        SIFT3D_IM_GET_UNITS(dst)[dim2] = SIFT3D_IM_GET_UNITS(src)[dim1];

	// Resize the output
	memcpy(SIFT3D_IM_GET_DIMS(dst), SIFT3D_IM_GET_DIMS(src), 
                IM_NDIMS * sizeof(int));
        SIFT3D_IM_GET_DIMS(dst)[dim1] = SIFT3D_IM_GET_DIMS(src)[dim2];
        SIFT3D_IM_GET_DIMS(dst)[dim2] = SIFT3D_IM_GET_DIMS(src)[dim1];
	dst->nc = src->nc;
	im_default_stride(dst);
	if (im_resize(dst))
		return SIFT3D_FAILURE;

	// Transpose the data
	SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)

                int src_coords[] = {x, y, z};
                int temp;

                // Permute the coordinates
                temp = src_coords[dim1];
                src_coords[dim1] = src_coords[dim2];
                src_coords[dim2] = temp;

                // Copy the datum
                SIFT3D_IM_GET_VOX(dst, x, y, z, c) = SIFT3D_IM_GET_VOX(src, 
                        src_coords[0], src_coords[1], src_coords[2], c);

	SIFT3D_IM_LOOP_END_C 

        return SIFT3D_SUCCESS;
}

/* Change an image's stride, preserving its data. 
 *
 * Parameters:
 *  -src: The source image. Must be initialized
 *  -strides: an array of length IM_NDIMS specifying the new strides
 *  -dst: The destination image. Must be initialized
 * 
 * Return: SIFT3D_SUCCESS (0) on success, nonzero otherwise. */
int im_restride(const Image * const src, const size_t *const strides,
		Image * const dst)
{

	int x, y, z, c;

	// Resize the output
	memcpy(SIFT3D_IM_GET_DIMS(dst), SIFT3D_IM_GET_DIMS(src), 
                IM_NDIMS * sizeof(int));
	memcpy(SIFT3D_IM_GET_STRIDES(dst), strides, IM_NDIMS * sizeof(size_t));
	dst->nc = src->nc;
	if (im_resize(dst))
		return SIFT3D_FAILURE;

	// Copy the data
	SIFT3D_IM_LOOP_START_C(dst, x, y, z, c)
	    SIFT3D_IM_GET_VOX(dst, x, y, z, c) =
	    SIFT3D_IM_GET_VOX(src, x, y, z, c);
	SIFT3D_IM_LOOP_END_C 

        return SIFT3D_SUCCESS;
}

/* Initializes a tform to ensure memory safety. 
 * Either this or the type-specific version must be called prior to using
 * a tform. */
int init_tform(void *const tform, const tform_type type)
{

	switch (type) {
	case TPS:
		puts("init_tform: TPS not yet implemented \n");
		return SIFT3D_FAILURE;
	case AFFINE:
		if (init_Affine((Affine *) tform, IM_NDIMS))
			return SIFT3D_FAILURE;
		break;
	default:
		puts("init_tform: unrecognized type \n");
		return SIFT3D_FAILURE;
	}

	return SIFT3D_SUCCESS;
}

/* Initialize an Affine struct. This initializes
 * all fields, and allocates memory for the inner
 * matrix, initializing it to zero. */
int init_Affine(Affine * const affine, const int dim)
{

	// Verify inputs
	if (dim < 2)
		return SIFT3D_FAILURE;

	// Initialize the type
	affine->tform.type = AFFINE;

	// Initialize the vtable
	affine->tform.vtable = &Affine_vtable;

	// Initialize the matrix
	if (init_Mat_rm(&affine->A, dim, dim + 1, SIFT3D_DOUBLE, SIFT3D_TRUE))
		return SIFT3D_FAILURE;

	return SIFT3D_SUCCESS;
}

/* Deep copy of a tform. Both src and dst must be initialized. */
int copy_tform(const void *const src, void *const dst)
{
	return TFORM_GET_VTABLE(src)->copy(src, dst);
}

/* Deep copy of one Affine to another. Both must be initialized. */
static int copy_Affine(const void *const src, void *const dst)
{

	const Affine *const srcAff = src;
	Affine *const dstAff = dst;

	return Affine_set_mat(&srcAff->A, dstAff);
}

/* Deep copy of one TPS to another. Both must be initialized. */
static int copy_Tps(const void *const src, void *const dst)
{
	SIFT3D_ERR("copy_Tps has not yet been implemented!");
	return SIFT3D_FAILURE;
}

/* Set an Affine transform to the given matrix.
 * mat is copied. mat must be an n x (n + 1) matrix, where
 * n is the dimensionality of the transformation. */
int Affine_set_mat(const Mat_rm * const mat, Affine * const affine)
{

	// Verify inputs
	if (mat->num_cols != mat->num_rows + 1 || mat->num_rows < 2)
		return SIFT3D_FAILURE;

	return convert_Mat_rm(mat, &affine->A, SIFT3D_DOUBLE);
}

/* Apply an arbitrary transformation to an [x, y, z] triple. */
void apply_tform_xyz(const void *const tform, const double x_in, 
                     const double y_in, const double z_in, double *const x_out,
		     double *const y_out, double *const z_out)
{
	TFORM_GET_VTABLE(tform)->apply_xyz(tform, x_in, y_in, z_in,
					   x_out, y_out, z_out);
}

/* Apply an Affine transformation to an [x, y, z] triple. */
static void apply_Affine_xyz(const void *const affine, const double x_in,
			     const double y_in, const double z_in,
			     double *const x_out, double *const y_out,
			     double *const z_out)
{

	const Affine *const aff = affine;

	const Mat_rm *const A = &aff->A;
	assert(AFFINE_GET_DIM(aff) == 3);
	*x_out = SIFT3D_MAT_RM_GET(A, 0, 0, double) * x_in +
	    SIFT3D_MAT_RM_GET(A, 0, 1, double) * y_in +
	    SIFT3D_MAT_RM_GET(A, 0, 2, double) * z_in +
	    SIFT3D_MAT_RM_GET(A, 0, 3, double);
	*y_out = SIFT3D_MAT_RM_GET(A, 1, 0, double) * x_in +
	    SIFT3D_MAT_RM_GET(A, 1, 1, double) * y_in +
	    SIFT3D_MAT_RM_GET(A, 1, 2, double) * z_in +
	    SIFT3D_MAT_RM_GET(A, 1, 3, double);
	*z_out = SIFT3D_MAT_RM_GET(A, 2, 0, double) * x_in +
	    SIFT3D_MAT_RM_GET(A, 2, 1, double) * y_in +
	    SIFT3D_MAT_RM_GET(A, 2, 2, double) * z_in +
	    SIFT3D_MAT_RM_GET(A, 2, 3, double);
}

/* Apply a thin-plate spline transformation to an [x, y, z] triple. */
static void apply_Tps_xyz(const void *const tps, const double x_in,
			  const double y_in, const double z_in,
			  double *const x_out, double *const y_out,
			  double *const z_out)
{

	const Tps *const t = tps;

	const Mat_rm *const params = &t->params;
	const Mat_rm *const kp_src = &t->kp_src;
	assert(t->dim == 3);
	int n;
	int ctrl_pts = kp_src->num_rows;	//number of control points
	double x_c, y_c, z_c, r_sq, U;
	double temp_x = 0.0, temp_y = 0.0, temp_z = 0.0;

	for (n = 0; n < ctrl_pts; n++) {
		x_c = SIFT3D_MAT_RM_GET(kp_src, n, 0, double);
		y_c = SIFT3D_MAT_RM_GET(kp_src, n, 1, double);
		z_c = SIFT3D_MAT_RM_GET(kp_src, n, 2, double);
		r_sq =
		    (x_in - x_c) * (x_in - x_c) + (y_in - y_c) * (y_in - y_c) +
		    (z_in - z_c) * (z_in - z_c);
		if (r_sq == 0) {
			U = 0.0;
		} else {
			U = r_sq * log(r_sq);
		}
		temp_x += U * SIFT3D_MAT_RM_GET(params, 0, n, double);
		temp_y += U * SIFT3D_MAT_RM_GET(params, 1, n, double);
		temp_z += U * SIFT3D_MAT_RM_GET(params, 2, n, double);
	}

	temp_x += SIFT3D_MAT_RM_GET(params, 0, ctrl_pts, double);
	temp_x += SIFT3D_MAT_RM_GET(params, 0, ctrl_pts + 1, double) * x_in;
	temp_x += SIFT3D_MAT_RM_GET(params, 0, ctrl_pts + 2, double) * y_in;
	temp_x += SIFT3D_MAT_RM_GET(params, 0, ctrl_pts + 3, double) * z_in;

	temp_y += SIFT3D_MAT_RM_GET(params, 1, ctrl_pts, double);
	temp_y += SIFT3D_MAT_RM_GET(params, 1, ctrl_pts + 1, double) * x_in;
	temp_y += SIFT3D_MAT_RM_GET(params, 1, ctrl_pts + 2, double) * y_in;
	temp_y += SIFT3D_MAT_RM_GET(params, 1, ctrl_pts + 3, double) * z_in;

	temp_z += SIFT3D_MAT_RM_GET(params, 2, ctrl_pts, double);
	temp_z += SIFT3D_MAT_RM_GET(params, 2, ctrl_pts + 1, double) * x_in;
	temp_z += SIFT3D_MAT_RM_GET(params, 2, ctrl_pts + 2, double) * y_in;
	temp_z += SIFT3D_MAT_RM_GET(params, 2, ctrl_pts + 3, double) * z_in;

	//Save results
	x_out[0] = temp_x;
	y_out[0] = temp_y;
	z_out[0] = temp_z;

}

/* Apply an arbitrary transform to a matrix. See apply_Affine_Mat_rm for
 * matrix formats. */
int apply_tform_Mat_rm(const void *const tform, const Mat_rm * const mat_in,
		       Mat_rm * const mat_out)
{
	return TFORM_GET_VTABLE(tform)->apply_Mat_rm(tform, mat_in, mat_out);
}

/* Apply a spline transformation to a matrix, 
 * by multiplication. See apply_Affine_Mat_rm for format of input matrices
 *
 * All matrices must be initialized with init_Mat_rm prior to use. For 3D!*/
static int apply_Tps_Mat_rm(const void *const tps, const Mat_rm * const mat_in,
			    Mat_rm * const mat_out)
{

	const Tps *const t = tps;

	//Spline transformation matrix is dim * [number of chosen points+dim+1]
	//sp_src is [number of chosen points] * dim
	const Mat_rm *const params = &(t->params);
	const Mat_rm *const kp_src = &(t->kp_src);

	int num_pts = mat_in->num_cols;	//number of points to be transformed
	int ctrl_pts = kp_src->num_rows;	//number of control points
	int m, n, q;
	double temp, x, y, z, r_sq, x_c, y_c, z_c;
	double U[ctrl_pts];
	//for each point
	for (q = 0; q < num_pts; q++) {
		//extract the coordinates
		x = SIFT3D_MAT_RM_GET(mat_in, 0, q, double);
		y = SIFT3D_MAT_RM_GET(mat_in, 1, q, double);
		z = SIFT3D_MAT_RM_GET(mat_in, 2, q, double);
		//Calculate U function for each control point
		for (n = 0; n < ctrl_pts; n++) {
			x_c = SIFT3D_MAT_RM_GET(kp_src, n, 0, double);
			y_c = SIFT3D_MAT_RM_GET(kp_src, n, 1, double);
			z_c = SIFT3D_MAT_RM_GET(kp_src, n, 2, double);
			r_sq =
			    (x - x_c) * (x - x_c) + (y - y_c) * (y - y_c) + (z -
									     z_c)
			    * (z - z_c);
			if (r_sq == 0) {
				U[n] = 0.0;
			} else {
				U[n] = r_sq * log(r_sq);
			}
		}
		//for each dimension
		for (m = 0; m < 3; m++) {	//For 3D!
			temp = 0.0;
			for (n = 0; n < ctrl_pts; n++) {
				temp +=
				    U[n] * SIFT3D_MAT_RM_GET(params, m, n,
							     double);
			}
			temp += SIFT3D_MAT_RM_GET(params, m, ctrl_pts, double);
			temp +=
			    SIFT3D_MAT_RM_GET(params, m, ctrl_pts + 1,
					      double) * x;
			temp +=
			    SIFT3D_MAT_RM_GET(params, m, ctrl_pts + 2,
					      double) * y;
			temp +=
			    SIFT3D_MAT_RM_GET(params, m, ctrl_pts + 3,
					      double) * z;

			//Store results
			SIFT3D_MAT_RM_GET(mat_out, m, q, double) = temp;
		}

	}
	return SIFT3D_SUCCESS;
}

/* Get the type of a tform. */
tform_type tform_get_type(const void *const tform)
{
	return ((Affine *) tform)->tform.type;
}

/* Get the size of a tform. */
size_t tform_get_size(const void *const tform)
{
	return TFORM_GET_VTABLE(tform)->get_size();
}

/* Get the size of a type of tform. */
size_t tform_type_get_size(const tform_type type)
{
	switch (type) {
	case AFFINE:
		return Affine_vtable.get_size();
	case TPS:
		return Tps_vtable.get_size();
	default:
		SIFT3D_ERR("tform_type_get_size: unrecognized " "type \n");
		return 0;
	}
}

/* Returns the size of an Affine struct */
static size_t Affine_get_size(void)
{
	return sizeof(Affine);
}

/* Returns the size of a Tps struct */
static size_t Tps_get_size(void)
{
	return sizeof(Tps);
}

/* Write a tform to a file. */
int write_tform(const char *path, const void *const tform)
{
	return TFORM_GET_VTABLE(tform)->write(path, tform);
}

/* Write an affine transformation to a file. */
static int write_Affine(const char *path, const void *const tform)
{

	const Affine *const affine = tform;

	return write_Mat_rm(path, &affine->A);
}

/* Write a thin-plate spline transformation to a file. */
static int write_Tps(const char *path, const void *const tform)
{
SIFT3D_IGNORE_UNUSED
	const Tps *const tps = tform;

	SIFT3D_ERR("write_Tps: this function has not yet been implemented.");
	return SIFT3D_FAILURE;
}

/* Free the memory associated with a tform */
void cleanup_tform(void *const tform)
{
	TFORM_GET_VTABLE(tform)->cleanup(tform);
}

/* Free the memory associated with an Affine transformation. */
static void cleanup_Affine(void *const affine)
{

	Affine *const aff = affine;

	cleanup_Mat_rm(&aff->A);
}

/* Free the memory assocaited with a thin-plate spline. */
static void cleanup_Tps(void *const tps)
{

	Tps *const t = tps;

	cleanup_Mat_rm(&t->params);
	cleanup_Mat_rm(&t->kp_src);
}

/* Apply an Affine transformation to a matrix, by multiplication. The format
 * of Mat_in should be:
 * [x1 x2 ... xN
 *  y1 y2 ... yN
 *      ...
 *  w1 w2 ... wN
 *  1  1  ... 1] 
 * 
 * mat_out will be resized to the appropriate size. The format will be:
 * [x1' x2' ... xN'
 *  y1' y2' ... yN'
 *      ...
 *  w1' w2' ... wN'] 
 *
 * All matrices must be initialized with init_Mat_rm prior to use. */
static int apply_Affine_Mat_rm(const void *const affine, 
        const Mat_rm * const mat_in, Mat_rm * const mat_out)
{

	const Affine *const aff = affine;

	return mul_Mat_rm(&aff->A, mat_in, mat_out);
}

/* Computes mat_in1 * mat_in2 = mat_out. mat_out will be resized
 * the appropriate size.
 *
 * All matrices must be initialized with init_Mat_rm prior to use. */
int mul_Mat_rm(const Mat_rm * const mat_in1, const Mat_rm * const mat_in2,
	       Mat_rm * const mat_out)
{

	int i, j, k;

	// Verify inputs
	if (mat_in1->num_cols != mat_in2->num_rows ||
	    mat_in1->type != mat_in2->type)
		return SIFT3D_FAILURE;

	// Resize mat_out
	mat_out->type = mat_in1->type;
	mat_out->num_rows = mat_in1->num_rows;
	mat_out->num_cols = mat_in2->num_cols;
	if (resize_Mat_rm(mat_out))
		return SIFT3D_FAILURE;

#define MAT_RM_MULTIPLY(type) \
    SIFT3D_MAT_RM_LOOP_START(mat_out, i, j) \
    type acc = 0; \
    for (k = 0; k < mat_in1->num_cols; k++) { \
    acc += SIFT3D_MAT_RM_GET(mat_in1, i, k, type) * \
    SIFT3D_MAT_RM_GET(mat_in2, k, j, type); \
    } \
    SIFT3D_MAT_RM_GET(mat_out, i, j, type) = acc; \
    SIFT3D_MAT_RM_LOOP_END

	// Row-major multiply
	switch (mat_out->type) {
	case SIFT3D_DOUBLE:
		MAT_RM_MULTIPLY(double) break;
	case SIFT3D_FLOAT:
		MAT_RM_MULTIPLY(float) break;
	case SIFT3D_INT:
		MAT_RM_MULTIPLY(int) break;
	default:
		puts("mul_Mat_rm: unknown type \n");
		return SIFT3D_FAILURE;
	}
#undef MAT_RM_MULTIPLY

	return SIFT3D_SUCCESS;
}

/* Computes the eigendecomposition of a real symmetric matrix, 
 * A = Q * diag(L) * Q', where Q is a real orthogonal matrix and L is a real 
 * diagonal matrix.
 *
 * A must be an [nxn] matrix. Q is [nxm], where m is in the interval [1, n],
 * depending on the values of A. L is [nx1], where the first m elements are
 * sorted in ascending order. The remaining n - m elements are zero. 
 * 
 * If Q is NULL, the eigenvectors will not be computed.
 *
 * The eigendecomposition is computed by divide and conquer.
 * 
 * This function resizes all non-null outputs and sets their type to double.
 *
 * This function does not ensure that A is symmetric.
 *
 * All matrices must be initialized prior to calling this funciton.
 * All matrices must have type double.
 *
 * Note: This function computes all of the eigenvalues, to a high degree of 
 * accuracy. A faster implementation is possible if you do not need high
 * precision, or if you do not need all of the eigenvalues, or if you do not 
 * need eigenvalues outside of some interval. 
 */
int eigen_Mat_rm(Mat_rm * A, Mat_rm * Q, Mat_rm * L)
{

	Mat_rm A_trans;
	double *work;
	fortran_int *iwork;
	double lwork_ret;
	fortran_int info, lwork, liwork;

	const char jobz = Q == NULL ? 'N' : 'V';
	const char uplo = 'U';
	const fortran_int n = A->num_cols;
	const fortran_int lda = n;
	const fortran_int lwork_query = -1;
	const fortran_int liwork_query = -1;

	// Verify inputs
	if (A->num_rows != n) {
		puts("eigen_Mat_rm: A be square \n");
		return SIFT3D_FAILURE;
	}
	if (A->type != SIFT3D_DOUBLE) {
		puts("eigen_Mat_rm: A must have type double \n");
		return SIFT3D_FAILURE;
	}
	// Resize outputs
	L->num_rows = n;
	L->num_cols = 1;
	L->type = SIFT3D_DOUBLE;
	if (resize_Mat_rm(L))
		return SIFT3D_FAILURE;

	// Initialize intermediate matrices and buffers
	work = NULL;
	iwork = NULL;
	if (init_Mat_rm(&A_trans, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE))
		goto EIGEN_MAT_RM_QUIT;

	// Copy the input matrix (A = A')
	if (copy_Mat_rm(A, &A_trans))
		goto EIGEN_MAT_RM_QUIT;

	// Query for the workspace sizes
	dsyevd_(&jobz, &uplo, &n, A_trans.u.data_double, &lda, L->u.data_double,
		&lwork_ret, &lwork_query, &liwork, &liwork_query, &info);

	if ((int32_t) info) {
		printf
		    ("eigen_Mat_rm: LAPACK dsyevd workspace query error code %d",
		     info);
		goto EIGEN_MAT_RM_QUIT;
	}
	// Allocate work spaces 
	lwork = (fortran_int) lwork_ret;
	if ((work = (double *)malloc(lwork * sizeof(double))) == NULL ||
	    (iwork =
	     (fortran_int *) malloc(liwork * sizeof(fortran_int))) == NULL)
		goto EIGEN_MAT_RM_QUIT;

	// Compute the eigendecomposition
	dsyevd_(&jobz, &uplo, &n, A_trans.u.data_double, &lda, L->u.data_double,
		work, &lwork, iwork, &liwork, &info);

	if ((int32_t) info) {
		printf("eigen_Mat_rm: LAPACK dsyevd error code %d", (int) info);
		goto EIGEN_MAT_RM_QUIT;
	}
	// Optionally return the eigenvectors
	if (Q != NULL && transpose_Mat_rm(&A_trans, Q))
		goto EIGEN_MAT_RM_QUIT;

	free(work);
	free(iwork);
	cleanup_Mat_rm(&A_trans);
	return SIFT3D_SUCCESS;

 EIGEN_MAT_RM_QUIT:
	if (work != NULL)
		free(work);
	if (iwork != NULL)
		free(iwork);
	cleanup_Mat_rm(&A_trans);
	return SIFT3D_FAILURE;
}

/* Solves the system AX=B exactly. A must be a square matrix.
 * This function first computes the reciprocal condition number of A.
 * If it is below the parameter "limit", it returns SIFT3D_SINGULAR. If limit 
 * is less than 0, a default value of 100 * eps is used.
 *
 * The system is solved by LU decomposition.
 * 
 * This function returns an error if A and B do not have valid dimensions. 
 * This function resizes X to [nx1] and changes the type to match B. 
 * All matrices must be initialized prior to calling this function.
 * All matrices must have type double.
 */
int solve_Mat_rm(const Mat_rm *const A, const Mat_rm *const B, 
        const double limit, Mat_rm *const X)
{

	Mat_rm A_trans, B_trans;
	double *work;
	fortran_int *ipiv, *iwork;
	double limit_arg, anorm, rcond;
	fortran_int info;

	const fortran_int m = A->num_rows;
	const fortran_int n = A->num_cols;
	const fortran_int nrhs = B->num_cols;
	const fortran_int lda = m;
	const fortran_int ldb = B->num_rows;
	const char norm_type = '1';
	const char trans = 'N';

	// Default parameters
	if (limit < 0)
		limit_arg = 100.0 * DBL_EPSILON;

	// Verify inputs
	if (m != n || ldb != m) {
		puts("solve_Mat_rm: invalid dimensions! \n");
		return SIFT3D_FAILURE;
	}
	if (A->type != SIFT3D_DOUBLE || B->type != SIFT3D_DOUBLE) {
		puts("solve_mat_rm: All matrices must have type double \n");
		return SIFT3D_FAILURE;
	}
	// Initialize intermediate matrices and buffers
	ipiv = NULL;
	work = NULL;
	iwork = NULL;
	if (init_Mat_rm(&A_trans, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE) ||
	    init_Mat_rm(&B_trans, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE) ||
	    (work = (double *)malloc(n * 4 * sizeof(double))) == NULL ||
	    (iwork = (fortran_int *) malloc(n * sizeof(fortran_int))) == NULL ||
	    (ipiv = (fortran_int *) calloc(m, sizeof(fortran_int))) == NULL)
		goto SOLVE_MAT_RM_QUIT;

	// Transpose matrices for LAPACK
	if (transpose_Mat_rm(A, &A_trans) || transpose_Mat_rm(B, &B_trans))
		goto SOLVE_MAT_RM_QUIT;

	// Compute the L1-norm of A
	anorm = dlange_(&norm_type, &m, &n, A_trans.u.data_double, &lda, work);

	// Compute the LU decomposition of A in place
	dgetrf_(&m, &n, A_trans.u.data_double, &lda, ipiv, &info);
	if ((int32_t) info < 0) {
		printf("solve_Mat_rm: LAPACK dgetrf error code %d \n", info);
		goto SOLVE_MAT_RM_QUIT;
	} else if ((int32_t) info > 0) {
		goto SOLVE_MAT_RM_SINGULAR;
	}
	// Compute the reciprocal condition number of A
	dgecon_(&norm_type, &n, A_trans.u.data_double, &lda, &anorm, &rcond,
		work, iwork, &info);
	if ((int32_t) info) {
		printf("solve_Mat_rm: LAPACK dgecon error code %d \n", info);
		goto SOLVE_MAT_RM_QUIT;
	}
	// Return if A is singular
	if (rcond < limit_arg)
		goto SOLVE_MAT_RM_SINGULAR;

	// Solve the system 
	dgetrs_(&trans, &n, &nrhs, A_trans.u.data_double, &lda, ipiv,
		B_trans.u.data_double, &ldb, &info);

	// Check for errors
	if ((int32_t) info) {
		printf("solve_Mat_rm: LAPACK dgetrs error code %d \n", info);
		goto SOLVE_MAT_RM_QUIT;
	}
	// Transpose results
	if (transpose_Mat_rm(&B_trans, X))
		goto SOLVE_MAT_RM_QUIT;

	free(ipiv);
	free(work);
	free(iwork);
	cleanup_Mat_rm(&A_trans);
	cleanup_Mat_rm(&B_trans);
	return SIFT3D_SUCCESS;

 SOLVE_MAT_RM_SINGULAR:
	free(ipiv);
	free(work);
	free(iwork);
	cleanup_Mat_rm(&A_trans);
	cleanup_Mat_rm(&B_trans);
	return SIFT3D_SINGULAR;

 SOLVE_MAT_RM_QUIT:
	if (ipiv != NULL)
		free(ipiv);
	if (work != NULL)
		free(work);
	if (iwork != NULL)
		free(iwork);
	cleanup_Mat_rm(&A_trans);
	cleanup_Mat_rm(&B_trans);
	return SIFT3D_FAILURE;
}

/* Solves the system AX=B by least-squares.
 *
 * A least-norm solution is computed using the singular 
 * value decomposition. A need not be full-rank.
 * 
 * This function returns an error if A and B do not have valid dimensions. 
 * This function resizes X to [nx1] and changes the type to match B. 
 * All matrices must be initialized prior to calling this funciton.
 * All matrices must have type double.
 */
int solve_Mat_rm_ls(const Mat_rm *const A, const Mat_rm *const B, 
        Mat_rm *const X)
{

	Mat_rm A_trans, B_trans;
	double *s, *work;
	double lwork_ret;
	fortran_int info, rank, lwork;
	int i, j;

	const double rcond = -1;
	const fortran_int m = A->num_rows;
	const fortran_int n = A->num_cols;
	const fortran_int nrhs = B->num_cols;
	const fortran_int lda = m;
	const fortran_int ldb = B->num_rows;
	const fortran_int lwork_query = -1;

	// Verify inputs 
	if (m != ldb) {
		puts("solve_Mat_rm_ls: invalid dimensions \n");
		return SIFT3D_FAILURE;
	}
	if (A->type != SIFT3D_DOUBLE || B->type != SIFT3D_DOUBLE) {
		puts("solve_mat_rm_ls: All matrices must have type double \n");
		return SIFT3D_FAILURE;
	}
	// Resize the output 
	X->type = SIFT3D_DOUBLE;
	X->num_rows = A->num_cols;
	X->num_cols = B->num_cols;
	if (resize_Mat_rm(X))
		return SIFT3D_FAILURE;

	// Initialize intermediate matrices and buffers
	s = NULL;
	work = NULL;
	if (init_Mat_rm(&A_trans, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE) ||
	    init_Mat_rm(&B_trans, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE) ||
	    (s = (double *)calloc(SIFT3D_MAX(m, n), sizeof(double))) == NULL)
		goto SOLVE_MAT_RM_LS_QUIT;

	// Transpose matrices for LAPACK
	if (transpose_Mat_rm(A, &A_trans) || transpose_Mat_rm(B, &B_trans))
		goto SOLVE_MAT_RM_LS_QUIT;

	// Get the size of the workspace
	dgelss_(&m, &n, &nrhs, A_trans.u.data_double, &lda,
		B_trans.u.data_double, &ldb, s, &rcond, &rank, &lwork_ret,
		&lwork_query, &info);
	if ((int32_t) info) {
		printf
		    ("solve_mat_rm: LAPACK dgelss work query error code %d \n",
		     info);
	}
	lwork = (fortran_int) lwork_ret;

	// Allocate the workspace
	if ((work = (double *)malloc(lwork * sizeof(double))) == NULL)
		goto SOLVE_MAT_RM_LS_QUIT;

	// Solve the system
	dgelss_(&m, &n, &nrhs, A_trans.u.data_double, &lda,
		B_trans.u.data_double, &ldb, s, &rcond, &rank, work, &lwork,
		&info);
	if ((int32_t) info) {
		printf("solve_mat_rm: LAPACK dgelss error code %d \n", info);
		goto SOLVE_MAT_RM_LS_QUIT;
	}
	// Transpose results to the new leading dimension
	SIFT3D_MAT_RM_LOOP_START(X, i, j)
	    SIFT3D_MAT_RM_GET(X, i, j, double) =
	    SIFT3D_MAT_RM_GET(&B_trans, j, i, double);
	SIFT3D_MAT_RM_LOOP_END free(s);
	free(work);
	cleanup_Mat_rm(&A_trans);
	cleanup_Mat_rm(&B_trans);
	return SIFT3D_SUCCESS;

 SOLVE_MAT_RM_LS_QUIT:
	if (s != NULL)
		free(s);
	if (work != NULL)
		free(work);
	cleanup_Mat_rm(&A_trans);
	cleanup_Mat_rm(&B_trans);
	return SIFT3D_FAILURE;
}

/* Computes the trace of a matrix. trace is assumed to be the same type as
 * mat. Returns an error if mat is not square. 
 * 
 * All matrices must be initialized with init_Mat_rm prior to calling 
 * this function. */
int trace_Mat_rm(Mat_rm * mat, void *trace)
{

	int i;

	// Verify inputs
	if (mat->num_rows != mat->num_cols || mat->num_rows < 1) {
		return SIFT3D_FAILURE;
	}
#define TRACE_MAT_RM(type) \
    {\
    type acc = 0; \
    for (i = 0; i < mat->num_rows; i++) { \
        acc += SIFT3D_MAT_RM_GET(mat, i, i, type); \
    } \
    *((type *) trace) = acc; \
    }

	// Take the trace
	switch (mat->type) {
	case SIFT3D_DOUBLE:
		TRACE_MAT_RM(double) break;
	case SIFT3D_FLOAT:
		TRACE_MAT_RM(float) break;
	case SIFT3D_INT:
		TRACE_MAT_RM(int) break;
	default:
		puts("trace_Mat_rm: unknown type \n");
		return SIFT3D_FAILURE;
	}
#undef TRACE_MAT_RM

	return SIFT3D_SUCCESS;
}

/* Tranposes a matrix. Resizes dst with the type of src. 
 * All matrices must be initialized prior to calling this function. */
int transpose_Mat_rm(const Mat_rm *const src, Mat_rm *const dst)
{

	int i, j;

	// Verify inputs
	if (src->num_rows < 1 || src->num_cols < 1)
		return SIFT3D_FAILURE;

	// Resize the output
	dst->type = src->type;
	dst->num_rows = src->num_cols;
	dst->num_cols = src->num_rows;
	if (resize_Mat_rm(dst))
		return SIFT3D_FAILURE;

#define TRANSPOSE_MAT_RM(type) \
    SIFT3D_MAT_RM_LOOP_START(src, i, j) \
    SIFT3D_MAT_RM_GET(dst, j, i, type) = \
                SIFT3D_MAT_RM_GET(src, i, j, type); \
    SIFT3D_MAT_RM_LOOP_END

	// Transpose
	switch (src->type) {
	case SIFT3D_DOUBLE:
		TRANSPOSE_MAT_RM(double);
		break;
	case SIFT3D_FLOAT:
		TRANSPOSE_MAT_RM(float);
		break;
	case SIFT3D_INT:
		TRANSPOSE_MAT_RM(int);
		break;
	default:
#ifndef NDEBUG
		puts("transpose_Mat_rm: unknown type \n");
#endif
		return SIFT3D_FAILURE;
	}
#undef TRANSPOSE_MAT_RM

	return SIFT3D_SUCCESS;
}

/* Computes the determinant of a symmetric matrix. det is assumed to be the 
 * same type as mat. Returns an error if mat is not square. 
 * 
 * This function does not verify that mat is symmetric.
 * 
 * All matrices must be initialized with init_Mat_rm prior to calling 
 * this function. */
int det_symm_Mat_rm(Mat_rm * mat, void *det)
{

	Mat_rm matd, L;
	double detd;
	int i, j;

	const int n = mat->num_cols;

	// Verify inputs
	if (n < 1 || mat->num_rows != n) {
		puts("det_symm_Mat_rm: invalid dimensions \n");
		return SIFT3D_FAILURE;
	}
	// Initialize intermediates
	if (init_Mat_rm(&matd, 0, 0, mat->type, SIFT3D_FALSE) ||
	    init_Mat_rm(&L, n, 1, SIFT3D_DOUBLE, SIFT3D_FALSE))
		goto DET_SYMM_QUIT;

	// Convert the matrix to type double
	if (convert_Mat_rm(mat, &matd, SIFT3D_DOUBLE))
		goto DET_SYMM_QUIT;

	// Get the eigendecomposition with LAPACK
	if (eigen_Mat_rm(&matd, NULL, &L))
		goto DET_SYMM_QUIT;

	// Take the determinant
	detd = 0.0;
	SIFT3D_MAT_RM_LOOP_START(&L, i, j)
	    detd += SIFT3D_MAT_RM_GET(&L, i, j, double);
	SIFT3D_MAT_RM_LOOP_END
	    // Convert the output to the correct type
	    switch (mat->type) {
	case SIFT3D_DOUBLE:
		*((double *)det) = detd;
		break;
	case SIFT3D_FLOAT:
		*((float *)det) = (float)detd;
		break;
	case SIFT3D_INT:
		*((int *)det) = (int)detd;
		break;
	default:
		puts("det_symm_Mat_rm: unknown type \n");
		goto DET_SYMM_QUIT;
	}

	cleanup_Mat_rm(&matd);
	cleanup_Mat_rm(&L);
	return SIFT3D_SUCCESS;

 DET_SYMM_QUIT:
	cleanup_Mat_rm(&matd);
	cleanup_Mat_rm(&L);
	return SIFT3D_FAILURE;
}

/* Apply a separable filter in multiple dimensions. This function resamples the
 * input to have the same units as f, then resamples the output to the
 * original units.
 *
 * Parameters:
 *  -src: The input image.
 *  -dst: The filtered image.
 *  -f: The filter to apply.
 *  -unit: The physical units of the filter kernel. Use -1.0 for the default,
 *      which is the same units as src.
 *
 * Return: SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int apply_Sep_FIR_filter(const Image * const src, Image * const dst,
			 Sep_FIR_filter * const f, const double unit)
{

	Image temp;
	Image *cur_src, *cur_dst;
	int i;

        const double unit_default = -1.0;

        // Verify inputs
        if (unit < 0 && unit != unit_default) {
                SIFT3D_ERR("apply_Sep_FIR_filter: invalid unit: %f, use "
                        "%f for default \n", unit, unit_default);
                return SIFT3D_FAILURE;
        }

        // Resize the output
        if (im_copy_dims(src, dst))
                return SIFT3D_FAILURE; 

	// Allocate temporary storage
	init_im(&temp);
	if (im_copy_data(src, &temp))
		goto apply_sep_f_quit;

#define SWAP_BUFFERS \
    if (cur_dst == &temp) { \
    cur_src = &temp; \
    cur_dst = dst; \
    } else { \
    cur_src = dst; \
    cur_dst = &temp; \
    }

	// Apply in n dimensions
	cur_src = (Image *) src;
	cur_dst = &temp;
	for (i = 0; i < IM_NDIMS; i++) {

                // Check for default parameters
                const double unit_arg = unit == unit_default ?
                        SIFT3D_IM_GET_UNITS(src)[i] : unit;

#ifdef SIFT3D_USE_OPENCL
                convolve_sep(cur_src, cur_dst, f, i, unit_arg);
		SWAP_BUFFERS
#else
                // Transpose so that the filter dimension is x
                if (i != 0) {
                        if (im_permute(cur_src, 0, i, cur_dst))
				goto apply_sep_f_quit;
                        SWAP_BUFFERS
                }

		// Apply the filter
		convolve_sep(cur_src, cur_dst, f, 0, unit_arg);
		SWAP_BUFFERS

                // Transpose back
                if (i != 0) {
			if (im_permute(cur_src, 0, i, cur_dst))
				goto apply_sep_f_quit;
                        SWAP_BUFFERS

		}
#endif
	}

	// Swap back
	SWAP_BUFFERS;

#undef SWAP_BUFFERS

	// Copy result to dst, if necessary
	if (cur_dst != dst && im_copy_data(cur_dst, dst))
		goto apply_sep_f_quit;

	// Clean up
	im_free(&temp);
	return SIFT3D_SUCCESS;

apply_sep_f_quit:
	im_free(&temp);
	return SIFT3D_FAILURE;
}

/* Initialize a separable FIR filter struct with the given parameters. If OpenCL
 * support is enabled and initialized, this creates a program to apply it with
 * separable filters.  
 *
 * Note that the kernel data will be copied, so the user can free it without 
 * affecting f. */
int init_Sep_FIR_filter(Sep_FIR_filter *const f, const int dim, const int width,
			const float *const kernel, const int symmetric)
{

        const size_t kernel_size = width * sizeof(float);

        // Save the data
	f->dim = dim;
	f->width = width;
	f->symmetric = symmetric;

        // Allocate the kernel memory
        if ((f->kernel = (float *) malloc(kernel_size)) == NULL) {
                SIFT3D_ERR("init_Sep_FIT_filter: out of memory! \n");
                return SIFT3D_FAILURE;
        }

        // Copy the kernel data
        memcpy(f->kernel, kernel, kernel_size);

#ifdef SIFT3D_USE_OPENCL
	{
		char src[1 << 15];
		char *template;
		cl_program program;
		cl_int err;
		float k;
		int i;

		const char *path = SEP_FIR_3D_PATH;
		const int half_width = f->half_width;

		// Load the template
		if ((template = read_file(path)) == NULL) {
			printf("init_Sep_FIR_Filter: error reading path %s \n",
			       path);
			return SIFT3D_FAILURE;
		}
		sprintf(src, "%s\n", template);

		// Write the unrolled kernel
		for (i = -half_width; i < half_width; i++) {
			k = f->kernel[i];
			sprintf(src, "acc += %.16f * "
				"read_imagef(src, sampler, center + d_xyz * %d); \n",
				k, i);
		}

		// Write the ending
		sprintf(src,
			"write_imagef(dst, sampler, (float4) center); \n } \n");

		// Compile the program  
		if (compile_cl_program_from_source(&program, cl_data.context,
						   cl_data.devices,
						   cl_data.num_devices,
						   (char **)&src, 1))
			return SIFT3D_FAILURE;
		f->cl_apply_unrolled =
		    clCreateKernel(program, "sep_fir_3d", &err);
		check_cl_error(err, "init_Sep_FIR_Filter: create kernel");
		clReleaseProgram(program);
	}
#endif
	return SIFT3D_SUCCESS;
}

/* Free a Sep_FIR_Filter. */
void cleanup_Sep_FIR_filter(Sep_FIR_filter *const f)
{

	if (f->kernel != NULL) {
		free(f->kernel);
		f->kernel = NULL;
	}
#ifdef SIFT3D_USE_OPENCL
	//TODO release OpenCL program
#endif
}

/* Initialize the values of im so that it can be used by the
 * resize function. Does not allocate memory. */
void init_im(Image *const im)
{
	im->data = NULL;
	im->cl_valid = SIFT3D_FALSE;

	im->ux = 1;
	im->uy = 1;
	im->uz = 1;

	im->size = 0;
	im->s = -1.0;
	memset(SIFT3D_IM_GET_DIMS(im), 0, IM_NDIMS * sizeof(int));
	memset(SIFT3D_IM_GET_STRIDES(im), 0, IM_NDIMS * sizeof(size_t));
}

/* Initialize a normalized Gaussian filter, of the given sigma.
 * If SIFT3D_GAUSS_WIDTH_FCTR is defined, use that value for
 * the ratio between the width of the filter and sigma. Otherwise,
 * use the default value 3.0 
 */
#ifndef SIFT3D_GAUSS_WIDTH_FCTR
#define SIFT3D_GAUSS_WIDTH_FCTR 3.0
#endif
int init_Gauss_filter(Gauss_filter * const gauss, const double sigma,
		      const int dim)
{

	float *kernel;
	double x;
	float acc;
	int i;

	const int half_width = sigma > 0 ? 
                SIFT3D_MAX((int)ceil(sigma * SIFT3D_GAUSS_WIDTH_FCTR), 1) :
                1;
	const int width = 2 * half_width + 1;

	// Initialize intermediates 
	if ((kernel = (float *) malloc(width * sizeof(float))) == NULL)
		return SIFT3D_FAILURE;

	// Calculate coefficients
	acc = 0;
	for (i = 0; i < width; i++) {
		// distance away from center of filter
		x = (double)i - half_width;

		// (x / sigma)^2 = x*x / (sigma*sigma)
		x /= sigma + DBL_EPSILON;

		// exponentiate result
		kernel[i] = (float)exp(-0.5 * x * x);

		// sum of all kernel elements
		acc += kernel[i];
	}

	// normalize kernel to sum to 1
	for (i = 0; i < width; i++) {
		kernel[i] /= acc;
	}

	// Save the filter data 
	gauss->sigma = sigma;
	if (init_Sep_FIR_filter(&gauss->f, dim, width, kernel,
				   SIFT3D_TRUE))
                goto init_Gauss_filter_quit;

        // Clean up
        free(kernel);

        return SIFT3D_SUCCESS;

init_Gauss_filter_quit:
        free(kernel);
        return SIFT3D_FAILURE;
}

/* Initialize a Gaussian filter to go from scale s_cur to s_next. */
int init_Gauss_incremental_filter(Gauss_filter * const gauss,
				  const double s_cur, const double s_next,
				  const int dim)
{
	double sigma;

	if (s_cur > s_next) {
                SIFT3D_ERR("init_Gauss_incremental_filter: "
                                "s_cur (%f) > s_next (%f) \n", s_cur, s_next);
                return SIFT3D_FAILURE;
        }
	assert(dim > 0);

	// Compute filter width parameter (sigma)
	sigma = sqrt(s_next * s_next - s_cur * s_cur);

	// Initialize filter kernel
	if (init_Gauss_filter(gauss, sigma, dim))
		return SIFT3D_FAILURE;

	return SIFT3D_SUCCESS;
}

/* Free a Gauss_filter */
void cleanup_Gauss_filter(Gauss_filter * gauss)
{
	cleanup_Sep_FIR_filter(&gauss->f);
}

/* Initialize a GSS filters stuct. This must be called before gss can be
 * used in any other functions. */
void init_GSS_filters(GSS_filters * const gss)
{
	gss->num_filters = -1;
	gss->gauss_octave = NULL;
}

/* Create GSS filters to create the given scale-space 
 * pyramid. */
int make_gss(GSS_filters * const gss, const Pyramid * const pyr)
{

	Image *cur, *next;
	int o, s;

	const int dim = 3;

	const int num_filters = pyr->num_levels - 1;
	const int first_level = pyr->first_level;
	const int last_level = SIFT3D_PYR_LAST_LEVEL(pyr);

	// Verify inputs
	if (num_filters < 1) {
		SIFT3D_ERR("make_gss: pyr has only %d levels, must have "
				"at least 2", pyr->num_levels);
		return SIFT3D_FAILURE;
	}

	// Free all previous data, if any
	cleanup_GSS_filters(gss);
	init_GSS_filters(gss);

	// Copy pyramid parameters
	gss->num_filters = num_filters;
	gss->first_level = first_level;

	// Allocate the filter array (num_filters cannot be zero)
	if ((gss->gauss_octave = (Gauss_filter *) 
		SIFT3D_safe_realloc(gss->gauss_octave, 
		num_filters * sizeof(Gauss_filter))) == NULL)
		return SIFT3D_FAILURE;

	// Make the filter for the very first blur
	next = SIFT3D_PYR_IM_GET(pyr, pyr->first_octave, first_level);
	if (init_Gauss_incremental_filter(&gss->first_gauss, pyr->sigma_n,
					  next->s, dim))
		return SIFT3D_FAILURE;

	// Make one octave of filters (num_levels - 1)
	o = pyr->first_octave;
	for (s = first_level; s < last_level; s++) {
		cur = SIFT3D_PYR_IM_GET(pyr, o, s);
		next = SIFT3D_PYR_IM_GET(pyr, o, s + 1);
		if (init_Gauss_incremental_filter(SIFT3D_GAUSS_GET(gss, s),
						  cur->s, next->s, dim))
			return SIFT3D_FAILURE;
	}

	return SIFT3D_SUCCESS;
}

/* Free all memory associated with the GSS filters. gss cannot be reused
 * unless it is reinitialized. */
void cleanup_GSS_filters(GSS_filters * const gss)
{

	int i;

	const int num_filters = gss->num_filters;

	// We are done if gss has no filters
	if (num_filters < 1)
		return;

	// Free the first filter
	cleanup_Gauss_filter(&gss->first_gauss);

	// Free the octave filters
	for (i = 0; i < num_filters; i++) {
		Gauss_filter *const g = gss->gauss_octave + i;
		cleanup_Gauss_filter(g);
	}

	// Free the octave filter buffer
	free(gss->gauss_octave);
}

/* Initialize a Pyramid for use. Must be called before a Pyramid can be used
 * in any other functions. */
void init_Pyramid(Pyramid * const pyr)
{
	pyr->levels = NULL;
        pyr->first_level = 0;
	pyr->num_levels = pyr->num_kp_levels = 0;
	pyr->first_octave = 0;
	pyr->num_octaves = 0;
        pyr->sigma0 = pyr->sigma_n = 0.0;
}

/* Resize a scale-space pyramid according to the size of base image im.
 *
 * Parameters:
 *  -im: An image with the desired dimensions and units at octave 0
 *  -first_level: The index of the first pyramid level per octave
 *  -num_kp_levels: The number of levels per octave in which keypoints are 
 *      detected
 *  -num_levels: The total number of levels. Must be greater than or equal to
 *      num_kp_levels.
 *  -first_octave: The index of the first octave (0 is the base)
 *  -num_octaves: The total number of octaves 
 *  -sigma0: The scale parameter of level 0, octave 0
 *  -sigma_n: The nominal scale of the image im.
 *  -pyr: The Pyramid to be resized.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int resize_Pyramid(const Image *const im, const int first_level, 
        const unsigned int num_kp_levels, const unsigned int num_levels,
        const int first_octave, const unsigned int num_octaves, 
        Pyramid *const pyr) {

        double units[IM_NDIMS];
        int dims[IM_NDIMS];
	double factor;
	int i, o, s;

	const double sigma0 = pyr->sigma0;
	const double sigma_n = pyr->sigma_n;
        const int old_num_total_levels = pyr->num_levels * pyr->num_octaves;
	const int num_total_levels = num_levels * num_octaves;

        // Verify inputs
        if (num_levels < num_kp_levels) {
                SIFT3D_ERR("resize_Pyramid: num_levels (%u) < "
                        "num_kp_levels (%d)", num_levels, num_kp_levels);
                return SIFT3D_FAILURE;
        }

        // Store the new parameters
        pyr->first_level = first_level;
        pyr->num_kp_levels = num_kp_levels;
        pyr->first_octave = first_octave;
        pyr->num_octaves = num_octaves;
        pyr->num_levels = num_levels;

        // Clean up old levels which are no longer needed 
        for (i = num_total_levels; i < old_num_total_levels; i++) {
                Image *const level = pyr->levels + i;
                im_free(level);
        }

	// Resize the outer array
        if (num_total_levels != 0 && 
		((pyr->levels = SIFT3D_safe_realloc(pyr->levels,
		num_total_levels * sizeof(Image))) == NULL))
                return SIFT3D_FAILURE;

	// We have nothing more to do if there are no levels
	if (num_total_levels == 0)
		return SIFT3D_SUCCESS;

        // Initalize new levels
        for (i = old_num_total_levels; i < num_total_levels; i++) {
                Image *const level = pyr->levels + i;
                init_im(level);
        }

        // We have nothing more to do if the image is empty
        if (im->data == NULL)
                return SIFT3D_SUCCESS;

	// Calculate base image dimensions and units
	factor = pow(2.0, -first_octave);
        for (i = 0; i < IM_NDIMS; i++) {
                dims[i] = (int) ((double) SIFT3D_IM_GET_DIMS(im)[i] * factor);
                units[i] = SIFT3D_IM_GET_UNITS(im)[i] * factor;
        }

	// Initialize each level separately
	SIFT3D_PYR_LOOP_START(pyr, o, s)
                        // Initialize Image fields
                        Image *const level = SIFT3D_PYR_IM_GET(pyr, o, s);
                        memcpy(SIFT3D_IM_GET_DIMS(level), dims, 
                                IM_NDIMS * sizeof(int));
                        memcpy(SIFT3D_IM_GET_UNITS(level), units, 
                                IM_NDIMS * sizeof(double));
	                level->nc = im->nc;
	                im_default_stride(level);

                        // Re-size data memory
                        if (im_resize(level))
                                return SIFT3D_FAILURE;

	        SIFT3D_PYR_LOOP_SCALE_END

	        // Adjust dimensions and recalculate image size
                for (i = 0; i < IM_NDIMS; i++) {
                        dims[i] /= 2;
                        units[i] *= 2;
                }

	SIFT3D_PYR_LOOP_OCTAVE_END 

        // Set the scales for the new levels
        return set_scales_Pyramid(pyr->sigma0, pyr->sigma_n, pyr);
}

/* Set the scale-space parameters on a Pyramid struct. Operates on all levels
 * of the pyramid. This function is called automatically by resize_Pyramid.
 *
 * Parameters:
 *  -sigma0: The scale parameter of level 0, octave 0
 *  -sigma_n: The nominal scale parameter of images being transfomed into
 *      this pyramid struct. 
 *  -Pyr: The Pyramid to be modified. */
int set_scales_Pyramid(const double sigma0, const double sigma_n, 
        Pyramid *const pyr) {

        int o, s;

        const int num_kp_levels = pyr->num_kp_levels;
        const Image *const first_level = 
                SIFT3D_PYR_IM_GET(pyr, pyr->first_octave, pyr->first_level);

        // Compute the scales of each level
        SIFT3D_PYR_LOOP_START(pyr, o, s)

                // Compute the scale 
                Image *const level = SIFT3D_PYR_IM_GET(pyr, o, s);
                const double scale = 
                        sigma0 * pow(2.0, o + (double) s / num_kp_levels);

                // Verify that sigma_n is not too large
                if (o == pyr->first_octave && s == pyr->first_level && 
                        scale < sigma_n) {
                        SIFT3D_ERR("set_scales_Pyramid: sigma_n too large "
                                "for these settings. Max allowed: %f \n", 
                                scale - DBL_EPSILON);
                        return SIFT3D_FAILURE;
                }

                // Save the scale
                level->s = scale;
        SIFT3D_PYR_LOOP_END

        // Store the parameters
        pyr->sigma0 = sigma0;
        pyr->sigma_n = sigma_n;

        return SIFT3D_SUCCESS;
}

/* Make a deep copy of a pyramid. */
int copy_Pyramid(const Pyramid * const src, Pyramid * const dst)
{

        Image dummy;
        const Image *base;
	int o, s, have_levels;

        // Initialize intermediates
        init_im(&dummy);

        // Set the scale parameters
        if (set_scales_Pyramid(src->sigma0, src->sigma_n, dst))
                return SIFT3D_FAILURE;

        // Get the base image
	if (src->levels == NULL || src->num_octaves <= 0 ||
	    src->num_levels <= 0) {
                base = &dummy;
                have_levels = SIFT3D_FALSE;
        } else {
                base = src->levels;
                have_levels = SIFT3D_TRUE;
        }

        // Resize dst
        if (resize_Pyramid(base, src->first_level, src->num_kp_levels,
                src->num_levels, src->first_octave, src->num_octaves, dst))
                goto copy_Pyramid_failure;

        // We are done if src has no levels
        if (!have_levels)
                goto copy_Pyramid_success;

	// Copy the levels
	SIFT3D_PYR_LOOP_START(dst, o, s)

		const Image *const src_level = SIFT3D_PYR_IM_GET(src, o, s);
		Image *const dst_level = SIFT3D_PYR_IM_GET(dst, o, s);

		if (src_level->data != NULL && 
                        im_copy_data(src_level, dst_level))
			return SIFT3D_FAILURE;

	SIFT3D_PYR_LOOP_END 

copy_Pyramid_success:
        im_free(&dummy);
	return SIFT3D_SUCCESS;

copy_Pyramid_failure:
        im_free(&dummy);
	return SIFT3D_FAILURE;
}

/* Release all memory associated with a Pyramid. pyr cannot be used again,
 * unless it is reinitialized. */
void cleanup_Pyramid(Pyramid * const pyr)
{

	int o, s;

	// We are done if there are no levels
	if (pyr->levels == NULL)
		return;

	// Free the levels
	SIFT3D_PYR_LOOP_START(pyr, o, s)
		Image *const level = SIFT3D_PYR_IM_GET(pyr, o, s);
		im_free(level);
	SIFT3D_PYR_LOOP_END

	// Free the pyramid level buffer
	free(pyr->levels);
}

/* Initialize a Slab for first use */
void init_Slab(Slab *const slab) {
    slab->buf_size = slab->num = 0;
    slab->buf = NULL;
}

/* Free all memory associated with a slab. Slab cannot be re-used after 
 * calling this function, unless re-initialized. */
void cleanup_Slab(Slab * const slab)
{
        if (slab->buf != NULL)
	        free(slab->buf);
}

/* Write the levels of a pyramid to separate files
 * for debugging. The path is prepended to the
 * octave and scale number of each image. 
 *
 * File type is inferred from the extension in path.
 *
 * Supported file formats:
 * -NIFTI
 */
int write_pyramid(const char *path, Pyramid * pyr)
{

	char path_appended[1024];
	int o, s;

	// Validate or create output directory
	if (mkpath(path, out_mode))
		return SIFT3D_FAILURE;

	// Save each image a separate file
	SIFT3D_PYR_LOOP_START(pyr, o, s)
	    sprintf(path_appended, "%s_o%i_s%i", path, o, s);
	if (write_nii(path_appended, SIFT3D_PYR_IM_GET(pyr, o, s)))
		return SIFT3D_FAILURE;
	SIFT3D_PYR_LOOP_END return SIFT3D_SUCCESS;
}

/* Exit and print a message to stdout. */
void err_exit(const char *str)
{
	SIFT3D_ERR("Error! Exiting at %s \n", str);
	exit(1);
}

/* Read a whole ASCII file into a string. Returns NULL
 * on error. */
SIFT3D_IGNORE_UNUSED
static char *read_file(const char *path)
{

	FILE *file;
	char *buf;
	size_t len;

	if ((file = fopen(path, "r")) == NULL) {
		return NULL;
	}
	fseek(file, 0, SEEK_END);
	len = ftell(file);
	rewind(file);
	if (ferror(file) || ((buf = malloc(len)) == NULL))
		return NULL;
	fread(buf, sizeof(char), len, file);

	return ferror(file) ? NULL : buf;
}

/* Ensure all directories in the given path exist.
 * Thanks to Jonathan Leffler
 * Modifications: Ignore everything after the last '/' 
 */
static int mkpath(const char *path, mode_t mode)
{
	char *pp, *sp, *copypath;
	int status;

	if ((copypath = strndup(path, FILENAME_MAX)) == NULL)
		status = -1;

	/* Ignore everything after the last '/' */
	if ((sp = strrchr(copypath, '/')) != NULL) {
		*sp = '\0';
	} else {
		/* If there is no '/', we have nothing to do */
		free(copypath);
		return SIFT3D_SUCCESS;
	}

	status = 0;
	pp = copypath;
	while (status == 0 && (sp = strchr(pp, '/')) != NULL) {
		if (sp != pp) {
			/* Neither root nor double slash in path */
			*sp = '\0';
			status = do_mkdir(copypath, mode);
			*sp = '/';
		}
		pp = sp + 1;
	}
	if (status == 0)
		status = do_mkdir(copypath, mode);

	free(copypath);
	return (status);
}

/* Make a directory if it does not exist.
 * Thanks to Jonathan Leffler */
static int do_mkdir(const char *path, mode_t mode)
{
	struct stat st;
	int status = 0;

	if (stat(path, &st) != 0) {
		/* Directory does not exist. EEXIST for race condition */
		if (cross_mkdir(path, mode) != 0 && errno != EEXIST)
			status = -1;
	} else if (!S_ISDIR(st.st_mode)) {
		errno = ENOTDIR;
		status = -1;
	}

	return (status);
}

/* Cross-platform mkdir */
static int cross_mkdir(const char *path, mode_t mode) {
#ifdef _MINGW_WINDOWS
        return mkdir(path);
#elif defined( _WINDOWS )
        return _mkdir(path);
#else
        return mkdir(path, mode);
#endif
}

/* Initialize a Tps struct. This initializes
 * all fields, and allocates memory for the inner
 * matrix, initializing it to zero. */
int init_Tps(Tps * tps, int dim, int terms)
{
	// Verify inputs
	if (dim < 2)
		return SIFT3D_FAILURE;

	// Initialize the type
	tps->tform.type = TPS;

	// Initialize the vtable
	tps->tform.vtable = &Tps_vtable;

	// Initialize the matrices
	if (init_Mat_rm(&tps->params, dim, terms, SIFT3D_DOUBLE, SIFT3D_TRUE))
		return SIFT3D_FAILURE;

	if (init_Mat_rm(&tps->kp_src, terms - dim - 1, dim,
			SIFT3D_DOUBLE, SIFT3D_TRUE))
		return SIFT3D_FAILURE;

	tps->dim = dim;
	return SIFT3D_SUCCESS;
}

/* Initialize a RANSAC struct with the default parameters */
void init_Ransac(Ransac *const ran)
{
	ran->err_thresh = SIFT3D_err_thresh_default;
	ran->num_iter = SIFT3D_num_iter_default;
}

/* Set the err_thresh parameter in a Ransac struct, checking for validity. */
int set_err_thresh_Ransac(Ransac *const ran, double err_thresh)
{

	if (err_thresh < 0.0) {
		SIFT3D_ERR("set_err_thresh_Ransac: invalid error "
			"threshold: %f \n", err_thresh);
		return SIFT3D_FAILURE;
	}

	ran->err_thresh = err_thresh;

	return SIFT3D_SUCCESS;
}

/* Set the num_iter parameter in a Ransac struct. */
int set_num_iter_Ransac(Ransac *const ran, int num_iter)
{
        if (num_iter < 1) {
                SIFT3D_ERR("set_num_iter_Ransac: invalid number of "
                                "iterations: %d \n", num_iter);
                return SIFT3D_FAILURE;
        }

	ran->num_iter = num_iter;

        return SIFT3D_SUCCESS;
}

/* Copy a Ransac struct from src to dst. */
int copy_Ransac(const Ransac *const src, Ransac *const dst) {
        return set_num_iter_Ransac(dst, src->num_iter) ||
                set_err_thresh_Ransac(dst, src->err_thresh);
}

/* Returns an array of k integers, (uniformly) randomly chosen from the 
 * integers 0 through n - 1.
 *
 * The value of *ret must either be NULL, or a pointer to a previously
 * allocated block. On successful return, *ret contains the k random integers.
 *
 * Returns SIFT3D_SUCCESS on succes, SIFT3D_FAILURE otherwise. */
static int n_choose_k(const int n, const int k, int **ret) {

        int i;

        // Verify inputs
        if (n < k || k < 1)
                goto n_choose_k_fail;

        // Allocate the array of n elements
        if ((*ret = malloc(n * sizeof(int))) == NULL)
                goto n_choose_k_fail;

        // Initialize the array of indices
        for (i = 0; i < n; i++) {
                (*ret)[i] = i;
        }

        // Randomize the first k indices using Knuth shuffles
        for (i = 0; i < k; i++) {

                int *const ints = *ret;
                const int temp = ints[i]; 
                const int rand_idx = i + rand() % (n - i);

                ints[i] = ints[rand_idx];
                ints[rand_idx] = temp;
        }

        // Release unused memory
        if ((*ret = SIFT3D_safe_realloc(*ret, k * sizeof(int))) == NULL)
                goto n_choose_k_fail;

        return SIFT3D_SUCCESS;

n_choose_k_fail:
        if (*ret != NULL) {
                free(*ret);
                *ret = NULL;
        }
        return SIFT3D_FAILURE;
}

//make the system matrix for spline
SIFT3D_IGNORE_UNUSED
static int make_spline_matrix(Mat_rm * src, Mat_rm * src_in, Mat_rm * sp_src,
			      int K_terms, int *r, int dim)
{
	int i, d;
	double x, y, z, x2, y2, z2, r_sq, U;
	src_in->type = SIFT3D_DOUBLE;
	sp_src->type = SIFT3D_DOUBLE;
	if (init_Mat_rm
	    (src_in, K_terms + dim + 1, K_terms + dim + 1, SIFT3D_DOUBLE,
	     SIFT3D_TRUE)) {
		return SIFT3D_FAILURE;
	}
	if (init_Mat_rm(sp_src, K_terms, dim, SIFT3D_DOUBLE, SIFT3D_TRUE)) {
		return SIFT3D_FAILURE;
	}
	for (i = 0; i < K_terms; i++) {
		//get the coordinate of current point
		switch (dim) {
		case 2:
			x = SIFT3D_MAT_RM_GET(src, r[i], 0, double);
			y = SIFT3D_MAT_RM_GET(src, r[i], 1, double);
			break;
		case 3:
			x = SIFT3D_MAT_RM_GET(src, r[i], 0, double);
			y = SIFT3D_MAT_RM_GET(src, r[i], 1, double);
			z = SIFT3D_MAT_RM_GET(src, r[i], 2, double);
			break;
		}
		for (d = 0; d < i; d++) {
			//compute r
			switch (dim) {
			case 2:
				x2 = SIFT3D_MAT_RM_GET(src, r[d], 0, double);
				y2 = SIFT3D_MAT_RM_GET(src, r[d], 1, double);
				r_sq =
				    (x - x2) * (x - x2) + (y - y2) * (y - y2);
				break;
			case 3:
				x2 = SIFT3D_MAT_RM_GET(src, r[d], 0, double);
				y2 = SIFT3D_MAT_RM_GET(src, r[d], 1, double);
				z2 = SIFT3D_MAT_RM_GET(src, r[d], 2, double);
				r_sq =
				    (x - x2) * (x - x2) + (y - y2) * (y - y2) +
				    (z - z2) * (z - z2);
				break;
			}
			//compute U
			U = r_sq * log(r_sq);
			//construct K
			SIFT3D_MAT_RM_GET(src_in, i, d, double) = U;
			SIFT3D_MAT_RM_GET(src_in, d, i, double) = U;
		}
		SIFT3D_MAT_RM_GET(src_in, i, i, double) = 0.0;
		//construct P and P'
		SIFT3D_MAT_RM_GET(src_in, i, K_terms, double) = 1.0;
		SIFT3D_MAT_RM_GET(src_in, K_terms, i, double) = 1.0;
		switch (dim) {
		case 2:
			SIFT3D_MAT_RM_GET(src_in, i, K_terms + 1, double) = x;
			SIFT3D_MAT_RM_GET(src_in, i, K_terms + 2, double) = y;
			SIFT3D_MAT_RM_GET(src_in, K_terms + 1, i, double) = x;
			SIFT3D_MAT_RM_GET(src_in, K_terms + 2, i, double) = y;
			break;
		case 3:
			SIFT3D_MAT_RM_GET(src_in, i, K_terms + 1, double) = x;
			SIFT3D_MAT_RM_GET(src_in, i, K_terms + 2, double) = y;
			SIFT3D_MAT_RM_GET(src_in, i, K_terms + 3, double) = z;
			SIFT3D_MAT_RM_GET(src_in, K_terms + 1, i, double) = x;
			SIFT3D_MAT_RM_GET(src_in, K_terms + 2, i, double) = y;
			SIFT3D_MAT_RM_GET(src_in, K_terms + 3, i, double) = z;
			break;
		}

		//construct sp_src matrix(matrix that stores control points)
		switch (dim) {
		case 2:
			SIFT3D_MAT_RM_GET(sp_src, i, 0, double) = x;
			SIFT3D_MAT_RM_GET(sp_src, i, 1, double) = y;
			break;
		case 3:
			SIFT3D_MAT_RM_GET(sp_src, i, 0, double) = x;
			SIFT3D_MAT_RM_GET(sp_src, i, 1, double) = y;
			SIFT3D_MAT_RM_GET(sp_src, i, 2, double) = z;
			break;
		}

	}

	//construct O
	for (i = 0; i < dim; i++) {
		for (d = 0; d < dim; d++) {
			SIFT3D_MAT_RM_GET(src_in, K_terms + i, K_terms + d,
					  double) = 0.0;
		}
	}

	return SIFT3D_SUCCESS;
}

//make the system matrix for affine
static int make_affine_matrix(const Mat_rm *const pts_in, const int dim, 
        Mat_rm *const mat_out)
{

	int i, j;

	const int num_rows = pts_in->num_rows;

	mat_out->type = SIFT3D_DOUBLE;
	mat_out->num_rows = num_rows;
	mat_out->num_cols = dim + 1;
	if (resize_Mat_rm(mat_out))
		return SIFT3D_FAILURE;

	for (i = 0; i < num_rows; i++) {

		//Add one row to the matrix
		for (j = 0; j < dim; j++) {
			SIFT3D_MAT_RM_GET(mat_out, i, j, double) =
			    SIFT3D_MAT_RM_GET(pts_in, i, j, double);
		}
		SIFT3D_MAT_RM_GET(mat_out, i, dim, double) = 1.0;
	}

	return SIFT3D_SUCCESS;
}

//extract the control matrix from tform struct (only valid for spline)
SIFT3D_IGNORE_UNUSED
static Mat_rm *extract_ctrl_pts(void *tform, tform_type type)
{
	Mat_rm *T;
	Tps *tps = (Tps *) tform;
	switch (type) {
	case TPS:
		T = extract_ctrl_pts_Tps(tps);
		break;
	case AFFINE:
		break;
	default:
		return NULL;
	}
	return T;
}

static Mat_rm *extract_ctrl_pts_Tps(Tps * tps)
{
	Mat_rm *kp_src = &tps->kp_src;
	return kp_src;
}

/* Solve for a transformation struct. 
 *
 * Paramters:
 *   src - See ransac().
 *   ref - See ransac()
 *   tform - See ransac()
 *
 * Returns SIFT3D_SUCCESS, SIFT3D_SINGULAR, or SIFT3D_FAILURE. See ransac() for
 * interpretation. */
static int solve_system(const Mat_rm *const src, const Mat_rm *const ref, 
        void *const tform)
{
	const tform_type type = tform_get_type(tform);

	//Mat_rm *kp_ref;
	Mat_rm ref_sys, X;
	int dim, ret;

	init_Mat_rm(&ref_sys, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE);
	init_Mat_rm(&X, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE);

	//construct source matrix and initialize reference vector
	switch (type) {
	case TPS:
		//kp_ref = extract_ctrl_pts(tform, type);
//          make_spline_matrix(ref, &ref_in, kp_ref, num_pts, r, dim);
		puts("solve_system: TPS not yet implemented");
		goto SOLVE_SYSTEM_FAIL;
	case AFFINE:
                dim = AFFINE_GET_DIM((Affine *const) tform);
		make_affine_matrix(ref, dim, &ref_sys);
		break;
	default:
		puts("solve_system: unknown type");
		goto SOLVE_SYSTEM_FAIL;
	}

	// solve for the coefficients                   
        ret = ref_sys.num_rows == ref_sys.num_cols ?
		solve_Mat_rm(&ref_sys, src, -1.0, &X) :
		solve_Mat_rm_ls(&ref_sys, src, &X);

	switch (ret) {
	case SIFT3D_SUCCESS:
		break;
	case SIFT3D_SINGULAR:
		goto SOLVE_SYSTEM_SINGULAR;
	default:
		goto SOLVE_SYSTEM_FAIL;
	}

	// Save the transformation matrix
	switch (type) {
	case TPS:
		//TODO
		goto SOLVE_SYSTEM_FAIL;
	case AFFINE:
        {
		Mat_rm X_trans;

		init_Mat_rm(&X_trans, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE);

		ret = transpose_Mat_rm(&X, &X_trans) ||
		    Affine_set_mat(&X_trans, (Affine *) tform);

		cleanup_Mat_rm(&X_trans);

		if (ret)
			goto SOLVE_SYSTEM_FAIL;

		break;
	}
	default:
		goto SOLVE_SYSTEM_FAIL;
	}

        // Clean up
	cleanup_Mat_rm(&ref_sys);
	cleanup_Mat_rm(&X);

	return SIFT3D_SUCCESS;

 SOLVE_SYSTEM_SINGULAR:
	cleanup_Mat_rm(&ref_sys);
	cleanup_Mat_rm(&X);
	return SIFT3D_SINGULAR;

 SOLVE_SYSTEM_FAIL:
	cleanup_Mat_rm(&ref_sys);
	cleanup_Mat_rm(&X);
	return SIFT3D_FAILURE;
}

//Find the SSD error for the i'th point
static double tform_err_sq(const void *const tform, const Mat_rm *const src, 
        const Mat_rm *const ref, const int i)
{

	double err = 0.0;
	//Initialization
	//in -- inputs coordinates of source points
	//out -- registered points
	//r -- reference points (ground truth)
	double x_in, y_in, z_in, x_r, y_r, z_r, x_out, y_out, z_out;

	//Find the source point
	x_in = SIFT3D_MAT_RM_GET(ref, i, 0, double);
	y_in = SIFT3D_MAT_RM_GET(ref, i, 1, double);
	z_in = SIFT3D_MAT_RM_GET(ref, i, 2, double);

	//Register
	apply_tform_xyz(tform, x_in, y_in, z_in, &x_out, &y_out, &z_out);

	//Find the reference point
	x_r = SIFT3D_MAT_RM_GET(src, i, 0, double);
	y_r = SIFT3D_MAT_RM_GET(src, i, 1, double);
	z_r = SIFT3D_MAT_RM_GET(src, i, 2, double);

	//Find the SSD error
	err = (x_r - x_out) * (x_r - x_out) + (y_r - y_out) * (y_r - y_out) +
	    (z_r - z_out) * (z_r - z_out);

	//return the result 
	return err;
}

/* Perform one iteration of RANSAC. 
 *
 * Parameters:
 *  src - The source points.
 *  ref - The reference points.
 *  tform - The output transformation. Must be initialized.
 *  cset - An array in which to store the concensus set. The value *cset must
 *         either be NULL, or a pointer to a previously allocated block.
 *  len - A location in which to store the length of the cset. 
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_SINGULAR if the system is 
 * near singular, and SIFT3D_FAILURE otherwise. */
static int ransac(const Mat_rm *const src, const Mat_rm *const ref, 
        const Ransac *const ran, void *tform, int **const cset, int *const len)
{
        int *rand_indices;
	Mat_rm src_rand, ref_rand;
	int i, j, num_rand, cset_len;

	const double err_thresh = ran->err_thresh;
	const double err_thresh_sq = err_thresh * err_thresh;
	const int num_pts = src->num_rows;
        const int num_dim = src->num_cols;
	const tform_type type = tform_get_type(tform);

	// Verify inputs
	if (src->type != SIFT3D_DOUBLE || src->type != ref->type) {
		puts("ransac: all matrices must have type double \n");
		return SIFT3D_FAILURE;
	}
	if (src->num_rows != ref->num_rows || src->num_cols != ref->num_cols) {
		puts("ransac: src and ref must have the same dimensions \n");
		return SIFT3D_FAILURE;
	}

        // Get the number of points for this transform
	switch (type) {
	case AFFINE:
		num_rand = AFFINE_GET_DIM((Affine *const) tform) + 1;
		break;
	default:
		printf("ransac: unknown transformation type \n");
                return SIFT3D_FAILURE;
	}

	// Initialize intermediates
        rand_indices = NULL;
	init_Mat_rm(&src_rand, num_rand, num_dim, SIFT3D_DOUBLE, SIFT3D_FALSE);
	init_Mat_rm(&ref_rand, num_rand, num_dim, SIFT3D_DOUBLE, SIFT3D_FALSE);

        // Draw random point indices
        if (n_choose_k(num_pts, num_rand, &rand_indices))
                goto RANSAC_FAIL;

        // Copy the random points
	SIFT3D_MAT_RM_LOOP_START(&src_rand, i, j)

                const int rand_idx = rand_indices[i];

                SIFT3D_MAT_RM_GET(&src_rand, i, j, double) =
                        SIFT3D_MAT_RM_GET(src, rand_idx, j, double);
                SIFT3D_MAT_RM_GET(&ref_rand, i, j, double) =
                        SIFT3D_MAT_RM_GET(ref, rand_idx, j, double);

        SIFT3D_MAT_RM_LOOP_END

        // Fit a transform to the random points
	switch (solve_system(&src_rand, &ref_rand, tform)) {
	case SIFT3D_SUCCESS:
		break;
	case SIFT3D_SINGULAR:
		goto RANSAC_SINGULAR;
	default:
		goto RANSAC_FAIL;
	}

	// Extract the consensus set
	cset_len = 0;
	for (i = 0; i < num_pts; i++) {

		// Calculate the error
		const double err_sq = tform_err_sq(tform, src, ref, i);

		// Reject points below the error threshold
		if (err_sq > err_thresh_sq)
			continue;

		// Add to the consensus set (++cset_len cannot be zero)
		if ((*cset = SIFT3D_safe_realloc(*cset, 
			++cset_len * sizeof(int))) == NULL)
			goto RANSAC_FAIL;

		(*cset)[cset_len - 1] = i;
	}

	// Return the new length of cset
	*len = cset_len;

        if (rand_indices != NULL)
                free(rand_indices);
	cleanup_Mat_rm(&src_rand);
	cleanup_Mat_rm(&ref_rand);
	return SIFT3D_SUCCESS;

RANSAC_SINGULAR:
        if (rand_indices != NULL)
                free(rand_indices);
	cleanup_Mat_rm(&src_rand);
	cleanup_Mat_rm(&ref_rand);
	return SIFT3D_SINGULAR;

RANSAC_FAIL:
        if (rand_indices != NULL)
                free(rand_indices);
	cleanup_Mat_rm(&src_rand);
	cleanup_Mat_rm(&ref_rand);
	return SIFT3D_FAILURE;
}

//Resize spline struct based on number of selected points
int resize_Tps(Tps * tps, int num_pts, int dim)
{
	Mat_rm *params = &(tps->params);
	Mat_rm *kp_src = &(tps->kp_src);
	params->num_cols = num_pts + dim + 1;
	params->num_rows = dim;
	kp_src->num_rows = num_pts;
	kp_src->num_cols = dim;
	if (resize_Mat_rm(params)) {
		return SIFT3D_FAILURE;
	}
	if (resize_Mat_rm(kp_src)) {
		return SIFT3D_FAILURE;
	}

	tps->dim = dim;
	return SIFT3D_SUCCESS;
}

/* Fit a transformation from ref to src points, using random sample concensus 
 * (RANSAC).
 * 
 * Parameters:
 *   ran - Struct storing RANSAC parameters.
 *   src - The [mxn] source points.
 *   ref - The [mxn] reference points.
 *   tform - The output transform. Must be initialized with init_from prior to 
 *           calling this function. 
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
int find_tform_ransac(const Ransac *const ran, const Mat_rm *const src, 
        const Mat_rm *const ref, void *const tform)
{

	Mat_rm ref_cset, src_cset;
	void *tform_cur;
	int *cset, *cset_best;
	int i, j, dim, num_terms, ret, len, len_best, min_num_inliers;

	const int num_iter = ran->num_iter;
	const int num_pts = src->num_rows;
	const size_t tform_size = tform_get_size(tform);
	const tform_type type = tform_get_type(tform);

	// Initialize data structures
	cset = cset_best = NULL;
	len_best = 0;
	if ((tform_cur = malloc(tform_size)) == NULL ||
	    init_tform(tform_cur, type) ||
	    init_Mat_rm(&src_cset, len_best, IM_NDIMS, SIFT3D_DOUBLE, 
            	SIFT3D_FALSE) ||
	    init_Mat_rm(&ref_cset, len_best, IM_NDIMS, SIFT3D_DOUBLE, 
		SIFT3D_FALSE))
		goto find_tform_quit;

	// initialize type-specific variables
	switch (type) {
	case AFFINE:
                dim = AFFINE_GET_DIM((Affine *const) tform);
		num_terms = dim + 1;
		min_num_inliers = 5;
		break;
	default:
		puts("find_tform_ransac: unsupported transformation "
		     "type \n");
		goto find_tform_quit;
	}

	if (num_pts < num_terms) {
		printf("Not enough matched points \n");
		goto find_tform_quit;
	}
	// Ransac iterations
	for (i = 0; i < num_iter; i++) {
		do {
			ret = ransac(src, ref, ran, tform_cur, &cset, &len);
		} while (ret == SIFT3D_SINGULAR);

		if (ret == SIFT3D_FAILURE)
			goto find_tform_quit;

		if (len > len_best) {
			len_best = len;
			if ((cset_best = (int *)SIFT3D_safe_realloc(cset_best,
				len * sizeof(int))) == NULL || 
				copy_tform(tform_cur, tform))
				goto find_tform_quit;
			memcpy(cset_best, cset, len * sizeof(int));
		}
	}

	// Check if the minimum number of inliers was found
	if (len_best < min_num_inliers) {
		puts("find_tform_ransac: No good model was found! \n");
		goto find_tform_quit; }

	// Resize the concensus set matrices
        src_cset.num_rows = ref_cset.num_rows = len_best;
        if (resize_Mat_rm(&src_cset) || resize_Mat_rm(&ref_cset))
                goto find_tform_quit;

	// Extract the concensus set
	SIFT3D_MAT_RM_LOOP_START(&src_cset, i, j)

	        const int idx = cset_best[i];

	        SIFT3D_MAT_RM_GET(&src_cset, i, j, double) =
	                SIFT3D_MAT_RM_GET(src, idx, j, double);
	        SIFT3D_MAT_RM_GET(&ref_cset, i, j, double) =
	                SIFT3D_MAT_RM_GET(ref, idx, j, double);

	SIFT3D_MAT_RM_LOOP_END

#ifdef SIFT3D_RANSAC_REFINE
	// Refine with least squares
	switch (solve_system(&src_cset, &ref_cset, tform_cur)) {
	case SIFT3D_SUCCESS:
		// Copy the refined transformation to the output
		if (copy_tform(tform_cur, tform))
			goto find_tform_quit;
		break;
	case SIFT3D_SINGULAR:
		// Stick with the old transformation 
#ifdef VERBOSE
		printf("find_tform_ransac: warning: least-squares refinement "
		       "abandoned due to numerical precision \n");
#endif
		break;
	default:
		goto find_tform_quit;
	}
#endif

        // Clean up
	free(cset);
	free(cset_best);
	cleanup_tform(tform_cur);
        cleanup_Mat_rm(&ref_cset);
        cleanup_Mat_rm(&src_cset);
	if (tform_cur != NULL)
		free(tform_cur);
	return SIFT3D_SUCCESS;

find_tform_quit:
        // Clean up and return an error
	if (cset != NULL)
		free(cset);
	if (cset_best != NULL)
		free(cset_best);
	cleanup_tform(tform_cur);
	if (tform_cur != NULL)
		free(tform_cur);
        cleanup_Mat_rm(&ref_cset);
        cleanup_Mat_rm(&src_cset);
	return SIFT3D_FAILURE;
}

/* Parse the GNU standard arguments (--version, --help). On return, the
 * getopt state is restored to the original.
 *
 * Return values:
 * -SIFT3D_HELP - "--help" was found
 * -SIFT3D_VERSION - "--version" was found, and the version message printed
 * -SIFT3D_FALSE - no GNU standard arguments were found */
int parse_gnu(const int argc, char *const *argv)
{

	int c;

	const int opterr_start = opterr;

	// Options
	const struct option longopts[] = {
		{"help", no_argument, NULL, SIFT3D_HELP},
		{"version", no_argument, NULL, SIFT3D_VERSION},
		{0, 0, 0, 0}
	};

	// Process the arguments
	opterr = 0;
	while ((c = getopt_long(argc, argv, "+", longopts, NULL)) != -1) {
		switch (c) {
		case SIFT3D_HELP:
			return SIFT3D_HELP;
		case SIFT3D_VERSION:
			puts(version_msg);
			return SIFT3D_VERSION;
		}
	}

	// Restore the state
	optind = 0;
	opterr = opterr_start;

	return SIFT3D_FALSE;
}

/* Print the bug message to stderr. */
void print_bug_msg()
{
	SIFT3D_ERR(bug_msg);
}
