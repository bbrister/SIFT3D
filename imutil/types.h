/* types.h
* ----------------------------------------------------------------
* Rice MRI Team
* ----------------------------------------------------------------
* This header contains data type definitions.
*-----------------------------------------------------------------
* Created: Blaine Rister 12/26/2013
* Last updated: Jenny Zhang 3/20/2014
*/

#include <time.h>

#ifndef _TYPES_H
#define _TYPES_H

#if _Win16 == 1 || _WIN32 == 1 || _WIN64 == 1 || \
	defined __WIN32__ || defined __TOS_WIN__ || \
	defined __WINDOWS__ && !defined _WINDOWS
#define _WINDOWS
#endif

/* OpenCL type definitions. */
#ifdef USE_OPENCL
#ifdef _WINDOWS
#include "CL\cl.h"
#else
#include "CL/cl.h"
#endif
#else
	/* Use dummy macro definitions */
#define CL_SUCCESS 0
#define CL_FAILURE -1

	/* Use dummy type definitions. */
	typedef unsigned char cl_uchar;
	typedef int cl_device_id;
	typedef int cl_command_queue;
	typedef int cl_device_type;
	typedef int cl_platform_id;
	typedef int cl_context;
	typedef int cl_image_format;
	typedef int cl_mem_flags;
	typedef int cl_kernel;
	typedef int cl_mem;
	typedef int cl_mem_object_type;
	typedef int cl_program;
	typedef int cl_int;
	typedef int cl_bool;
	typedef unsigned int cl_uint;
#endif

/* Parameters */
#define ORI_PER_KP 4	// Maximum number of orientations per keypoint
#define NBINS_AZ 8		// Number of bins for azimuthal angles
#define NBINS_PO 4		// Number of bins for polar angles
#define NHIST_PER_DIM 4 // Number of SIFT descriptor histograms per dimension 
#define ICOS_HIST			// Icosahedral gradient histogram

/* Constants */
#define IM_NDIMS 3 // Number of dimensions in an Image
#define ICOS_NFACES 20 // Number of faces in an icosahedron
#define ICOS_NVERT 12 // Number of vertices in an icosahedron

/* Derived constants */
#define DESC_NUM_TOTAL_HIST (NHIST_PER_DIM * NHIST_PER_DIM * NHIST_PER_DIM)
#define DESC_NUMEL (DESC_NUM_TOTAL_HIST * HIST_NUMEL)

// The number of elements in a gradient histogram
#ifdef ICOS_HIST
#define HIST_NUMEL (ICOS_NVERT)
#else
#define HIST_NUMEL (NBINS_AZ * NBINS_PO)
#endif

/* Possible data types for matrix elements */ 
typedef enum _data_type {
	DOUBLE,
	FLOAT,
	INT
} data_type;

/* Struct to hold OpenCL programs for this library */
typedef struct _kernels {
	cl_kernel downsample_2x_3d;
} Kernels;

/* Struct to hold OpenCL data about the user system */
typedef struct _CL_data {
	cl_device_id *devices;	  // num_devices elements
	cl_command_queue *queues; // One per device
	cl_platform_id platform;
	cl_context context;
	cl_uint num_devices;
	cl_image_format image_format;
	cl_mem_flags mem_flags;
	Kernels kernels;
	int valid;		// Is this struct valid?
} CL_data;

/* Struct to hold a dense matrix in row-major order */
typedef struct _Mat_rm {

	union {
		double *data_double;
		float  *data_float;
		int *data_int;
	} u;
	data_type type;
	int num_cols;
	int num_rows;	
	size_t numel;		// number of elements
	size_t size;		// size of the buffer, in bytes

} Mat_rm;

/* Struct to hold image data. The image is a rectangular prism, 
 * where the bottom-left corner is [0 0 0], the x-stride is 1,
 * the y-stride is the width in x, and the z-stride is the
 * size of an xy plane. For convenience use the macros IM_GET_IDX, 
 * IM_GET_VOX, and IM_SET_VOX to manipulate this struct. */
typedef struct _Image {

	float *data;		// Raster of voxel values ~16MB
	int *dims;		// Array-style access to {nx, ...}
	int *strides;		// Array-style access to {x_stride, ...}
	cl_mem cl_image;	// Same-sized OpenCL image object
	double s;		// scale-space location
	size_t size;		// Total size in pixels
	int nx, ny, nz;		// Dimensions in x, y, and z
	int x_stride;		// Stride in x direction
	int y_stride;		// Stride in y direction
	int z_stride;		// Stride in z direction
        int nc;                 // The number of channels
	int cl_valid;		// If TRUE, cl_image is valid

} Image;

/* Holds separable FIR filters and programs to apply them */
typedef struct _Sep_FIR_filter {

	cl_kernel cl_apply_unrolled;	// unrolled OpenCL program to apply filter
	float *kernel;	// filter weights
	int dim;		// dimensionality, e.g. 3 for MRI
	int width;		// number of weights				
	int half_width;	// floor(width / 2)
	int symmetric;	// enable symmetric optimizations: FALSE or TRUE

} Sep_FIR_filter;

/* Holds Gaussian filters */
typedef struct _Gauss_filter {

	double sigma;
	Sep_FIR_filter f;

} Gauss_filter;

/* Holds Gaussian Scale-Space filters */
typedef struct _GSS_filters {

	Gauss_filter first_gauss;		// Used on the very first blur
	Gauss_filter *gauss_octave;		// Array of kernels for one octave
	int num_filters;				// Number of filters for one octave
	int first_level;
	int last_level;

} GSS_filters;

/* Struct to hold SIFT detector filters */
typedef struct _SIFT_filters {

	GSS_filters gss;

} SIFT_filters;

/* Struct to hold miscellaneous SIFT detector OpenCL kernels */
typedef struct _SIFT_cl_kernels {

	cl_kernel downsample_2;

} SIFT_cl_kernels;

/* Struct to hold a scale-space image pyramid */
typedef struct _Pyramid {
	
	// Filters
	SIFT_filters filters;

	// Levels in all octaves
	Image *levels;	

	// Scale-space parameters
	double sigma_n;
	double sigma0;
	int num_kp_levels;

	// Indexing information -- see macros.h
	int num_octaves;
	int num_levels;
	int first_octave;
	int last_octave;
	int first_level;
	int last_level;	

} Pyramid;

/* Struct defining a vector in spherical coordinates */
typedef struct _Svec {

	float mag;	// Magnitude
	float po;	// Polar angle, [0, pi)
	float az;	// Azimuth angle, [0, 2pi)

} Svec;

/* Struct defining a vector in Cartesian coordinates */
typedef struct _Cvec {

	float x;
	float y;
	float z;

} Cvec;

/* Slab allocation struct */
typedef struct _Slab {

	void *buf;			// Buffer
	size_t num;			// Number of elements currently in buffer
	size_t buf_length;	// Maximum number of elements in buffer

} Slab;

/* Struct defining a keypoint in 3D space. */
typedef struct _Keypoint {

	float r_data[3 * 3];			// Memory for matrix R, do not use this
	Mat_rm R;						// Rotation matrix into Keypoint space
	double xd, yd, zd;				// sub-pixel x, y, z
	double  sd;						// absolute scale
	double  sd_rel;					// octave-relative scale
	int xi, yi, zi, o, s;			// [x, y, z], [o, s]

} Keypoint;

/* Struct to hold keypoints */
typedef struct _Keypoint_store {
	
	Keypoint *buf;
	Slab slab;
	int nx, ny, nz;		// dimensions of first octave

} Keypoint_store;

/* Struct defining an orientation histogram in
 * spherical coordinates. */
typedef struct _Hist {
	float bins[HIST_NUMEL];
} Hist;

/* Triangle */
typedef struct _Tri {
	Cvec v[3]; // Vertices
	int idx[3]; // Index of each vertex in the solid
} Tri;

/* Triangle mesh */
typedef struct _Mesh {
	Tri *tri; 	// Triangles
	int num;	// Number of triangles
} Mesh;

/* Struct defining a 3D SIFT descriptor */
typedef struct _SIFT3D_Descriptor {

	Hist hists[DESC_NUM_TOTAL_HIST]; // Array of orientation histograms
	double xd, yd, zd, sd;	// sub-pixel [x, y, z], absolute scale

} SIFT3D_Descriptor;

/* Struct to hold SIFT3D descriptors */
typedef struct _SIFT3D_Descriptor_store {

	SIFT3D_Descriptor *buf;
	size_t num;
	int nx, ny, nz;			// Image dimensions

} SIFT3D_Descriptor_store;

/* Struct to hold all parameters and internal data of the 
 * SIFT3D algorithms */
typedef struct _SIFT3D {

        // Triange mesh
	Mesh mesh;

	// Filters and filtering programs
	SIFT_filters filters;

	// Other OpenCL kernels
	SIFT_cl_kernels kernels;

	// Gaussian pyramid
	Pyramid gpyr;

	// DoG pyramid
	Pyramid dog;

	// Image to process
	Image *im;

	// Parameters
	double peak_thresh; // Keypoint peak threshold
	double corner_thresh; // Keypoint corner threshold
        double dense_sigma; // Gaussian window parameter for dense descriptors
        int dense_rotate; // If true, dense descriptors are rotation-invariant

	// Profiling info
	//timer_t start;

} SIFT3D;

/* Geometric transformations that can be applied by this library. */
typedef enum _tform_type {
	AFFINE,         // Affine (linear + constant)
	TPS             // Thin-plate spline	
} tform_type;

/* Interpolation algorithms that can be used by this library. */
typedef enum _interp_type {
        LINEAR,         // N-linear interpolation
        LANCZOS2        // Lanczos kernel, a = 2
} interp_type;

/* Struct to hold an affine transformation */
typedef struct _Affine {
	Mat_rm A;		// Transformation matrix, x' = Ax
	int dim; 		// Dimensionality, e.g. 3
} Affine;

/* Struct to hold a thin-plate spline */
typedef struct _TPS {
	Mat_rm params;	// Transformation matrix, dim * number of control point + dim +1
	Mat_rm kp_src;	// Control point matrix, number of control point * dim
	int dim; 	// Dimensionality, e.g. 3
} Tps;

/* Struct to hold RANSAC parameters */
typedef struct _Ransac {
	double min_inliers; //ratio of inliers threshold for RANSAC concensus sets
 	double err_thresh; //error threshold for RANSAC inliers
	int num_iter; //number of RANSAC iterations
} Ransac;

#endif
