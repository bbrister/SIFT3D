/* -----------------------------------------------------------------------------
 * imtypes.h
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This header contains data type definitions.
 * -----------------------------------------------------------------------------
 */

#include <time.h>

#ifndef _IMTYPES_H
#define _IMTYPES_H

#ifdef __cplusplus
extern "C" {
#endif

// Return codes
#define SIFT3D_SINGULAR 1
#define SIFT3D_SUCCESS 0
#define SIFT3D_FAILURE -1
#define SIFT3D_HELP 1
#define SIFT3D_VERSION 2

// Truth values
#define SIFT3D_TRUE 1
#define SIFT3D_FALSE 0

// Platform types
#if _Win16 == 1 || _WIN32 == 1 || _WIN64 == 1 || \
	defined __WIN32__ || defined __TOS_WIN__ || \
	defined __WINDOWS__ && !defined _WINDOWS
#define _WINDOWS
#endif
#if (defined(__MINGW32__) || defined(__MINGW64__)) && defined(_WINDOWS) && \
         !defined _MINGW_WINDOWS 
#define _MINGW_WINDOWS
#endif

/* Missing types for MEX */
#ifdef SIFT3D_MEX
#ifdef _MINGW_WINDOWS
// char16_t is not defined in MinGW with default settings
typedef wchar_t char16_t;
#endif
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

/* File separator character */
#ifdef _WINDOWS
#define SIFT3D_FILE_SEP '\\'
#else
#define SIFT3D_FILE_SEP '/'
#endif

/* Parameters */
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

/* Supported image file formats */
typedef enum _im_format {
        ANALYZE, /* Analyze */
        DICOM, /* DICOM */
        DIRECTORY, /* Directory */
        NIFTI, /* NIFTI-1 */ 
        UNKNOWN, /* Not one of the known extensions */
        ERROR /* Error occurred in determining the format */
} im_format;

/* Possible data types for matrix elements */ 
typedef enum _Mat_rm_type {
	DOUBLE,
	FLOAT,
	INT
} Mat_rm_type;

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
	size_t size;		// Size of the buffer, in bytes
	int num_cols;           // Number of columns 
	int num_rows;           // Number of rows	
        int static_mem;         // Flag for statically-allocated memory
	Mat_rm_type type;       // DOUBLE, FLOAT, or INT

} Mat_rm;

/* Struct to hold image data. The image is a rectangular prism, 
 * where the bottom-left corner is [0 0 0], the x-stride is 1,
 * the y-stride is the width in x, and the z-stride is the
 * size of an xy plane. For convenience use the macros IM_GET_IDX, 
 * IM_GET_VOX, and IM_SET_VOX to manipulate this struct. */
typedef struct _Image {

	float *data;		// Raster of voxel values ~16MB
	cl_mem cl_image;	// Same-sized OpenCL image object
	double s;		// scale-space location
	size_t size;		// Total size in pixels
	int nx, ny, nz;		// Dimensions in x, y, and z
	double ux, uy, uz;	// Real world dimensions in x, y, and z
        size_t xs, ys, zs;      // Stride in x, y, and z
        int nc;                 // The number of channels
	int cl_valid;		// If TRUE, cl_image is valid

} Image;

/* Holds separable FIR filters and programs to apply them */
typedef struct _Sep_FIR_filter {

	cl_kernel cl_apply_unrolled;	// unrolled OpenCL program to apply filter
	float *kernel;	// filter weights
	int dim;	// dimensionality, e.g. 3 for MRI
	int width;	// number of weights				
	int symmetric;	// enable symmetric optimizations: FALSE or TRUE

} Sep_FIR_filter;

/* Holds Gaussian filters */
typedef struct _Gauss_filter {

	double sigma;
	Sep_FIR_filter f;

} Gauss_filter;

/* Holds Gaussian Scale-Space filters */
typedef struct _GSS_filters {

	Gauss_filter first_gauss;	// Used on the very first blur
	Gauss_filter *gauss_octave;	// Array of kernels for one octave
	int num_filters;		// Number of filters for one octave
	int first_level;                // Index of the first scale level

} GSS_filters;

/* Struct to hold miscellaneous SIFT detector OpenCL kernels */
typedef struct _SIFT_cl_kernels {

	cl_kernel downsample_2;

} SIFT_cl_kernels;

/* Struct to hold a scale-space image pyramid */
typedef struct _Pyramid {
	
	// Levels in all octaves
	Image *levels;	

	// Scale-space parameters
	double sigma_n;
	double sigma0;
	int num_kp_levels;

	// Indexing information -- see immacros.h
	int first_octave;
	int num_octaves;
	int first_level;
	int num_levels;

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
	size_t buf_size;	        // Buffer capactiy, in bytes

} Slab;

/* Struct defining a keypoint in 3D space. */
typedef struct _Keypoint {

	float r_data[IM_NDIMS * IM_NDIMS];	// Memory for matrix R, do not use this
	Mat_rm R;				// Rotation matrix into Keypoint space
	double xd, yd, zd;			// sub-pixel x, y, z
	double  sd;				// absolute scale
	int o, s;			        // pyramid indices 

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

        // Filters for computing the GSS pyramid
	GSS_filters gss;

	// Other OpenCL kernels
	SIFT_cl_kernels kernels;

	// Gaussian pyramid
	Pyramid gpyr;

	// DoG pyramid
	Pyramid dog;

	// Image to process
	Image im;

	// Parameters
	double peak_thresh; // Keypoint peak threshold
	double corner_thresh; // Keypoint corner threshold
        int dense_rotate; // If true, dense descriptors are rotation-invariant

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

/* Virtual function table for Tform class */
typedef struct _Tform_vtable {

        int (*copy)(const void *const, void *const);

        void (*apply_xyz)(const void *const, const double, const double, 
                const double, double *const, double *const, double *const);

        int (*apply_Mat_rm)(const void *const, const Mat_rm *const, 
                Mat_rm *const);

        size_t (*get_size)(void);

        int (*write)(const char *, const void *const);
       
        void (*cleanup)(void *const);

} Tform_vtable;

/* "Abstract class" of transformations */
typedef struct _Tform {
        tform_type type; // The specific type, e.g. Affine, TPS
        const Tform_vtable *vtable; // Table of virtual functions
} Tform;

/* Struct to hold an affine transformation */
typedef struct _Affine {
        Tform tform;    // Abstract parent class
	Mat_rm A;	// Transformation matrix, x' = Ax
} Affine;

/* Struct to hold a thin-plate spline */
typedef struct _Tps {
        Tform tform;       // Abstract parent class
	Mat_rm params;	// Transformation matrix, dim * number of control point + dim +1
	Mat_rm kp_src;	// Control point matrix, number of control point * dim
	int dim; 	// Dimensionality, e.g. 3
} Tps;

/* Struct to hold RANSAC parameters */
typedef struct _Ransac {
 	double err_thresh; //error threshold for RANSAC inliers
	int num_iter; //number of RANSAC iterations
} Ransac;

#ifdef __cplusplus
}
#endif

#endif
