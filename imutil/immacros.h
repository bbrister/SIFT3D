/* -----------------------------------------------------------------------------
 * immacros.h
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This header defines preprocessor macros for the imutil library.
 * -----------------------------------------------------------------------------
 */

#include "imtypes.h"

#ifdef SIFT3D_MEX
#include <uchar.h>
#include "mex.h"
#else
#include "stdio.h"
#endif

#ifndef _IMMACROS_H
#define _IMMACROS_H

#ifdef __cplusplus
extern "C" {
#endif

// Print an error message
#ifdef SIFT3D_MEX
#define SIFT3D_ERR(...) \
        mexWarnMsgIdAndTxt("sift3d:internal", __VA_ARGS__)
#else
#define SIFT3D_ERR(...) fprintf(stderr, __VA_ARGS__)
#endif

// Math macros
#define SIFT3D_MIN(x, y) ((x) < (y) ? (x) : (y))
#define SIFT3D_MAX(x, y) ((x) > (y) ? (x) : (y))
#define SIFT3D_AZ_MAX_F (2 * (float) M_PI) // Maximum azimuth
#define SIFT3D_PO_MAX_F ((float) M_PI) // Maximum polar angle

// Compiler flags
#ifdef __GNUC__
#define SIFT3D_IGNORE_UNUSED __attribute__((unused))
#endif

// Get a pointer to the [nx, ny, nz] array of an image
#define SIFT3D_IM_GET_DIMS(im) \
        (&(im)->nx)

// Get a pointer to the [xs, ys, zs] array of an image
#define SIFT3D_IM_GET_STRIDES(im) \
        (&(im)->xs)

// Get a pointer to the [ux, uy, uz] array of an image
#define SIFT3D_IM_GET_UNITS(im) \
        (&(im)->ux)

// Get the index of an [x,y,z] pair in an image 
#define SIFT3D_IM_GET_IDX(im, x, y, z, c) ((size_t) (x) * (im)->xs + \
        (size_t) (y) * (im)->ys + (size_t) (z) * (im)->zs + (size_t) (c))

// Get the value of voxel [x,y,z] in an image 
#define SIFT3D_IM_GET_VOX(im, x, y, z, c) ((im)->data[ \
        SIFT3D_IM_GET_IDX((im), (x), (y), (z), (c))])

// Loop through an image in x, z, y order. Delmit with SIFT3D_IM_LOOP_END
#define SIFT3D_IM_LOOP_START(im, x, y, z) \
	for (z = 0; (z) < (im)->nz; (z)++) {	\
	for ((y) = 0; (y) < (im)->ny; (y)++) {	\
	for ((x) = 0; (x) < (im)->nx; (x)++) {

/* As in SIFT3D_IM_LOOP_START, but also loop through each channel */
#define SIFT3D_IM_LOOP_START_C(im, x, y, z, c) \
        SIFT3D_IM_LOOP_START(im, x, y, z) \
        for ((c) = 0; (c) < (im)->nc; (c)++) {

/* Loop through an image iterating with the (inclusive) x, y, z bounds given.
 * Delimit with SIFT3D_IM_LOOP_END. */
#define SIFT3D_IM_LOOP_LIMITED_START(im, x, y, z, x_start, x_end, \
			      y_start, y_end, z_start, z_end) \
	for (z = z_start; (z) <= z_end; (z)++) { \
	for ((y) = y_start; (y) <= y_end; (y)++) { \
	for ((x) = x_start; (x) <= x_end; (x)++) {		

/* As in SIFT3D_IM_LOOP_LIMITED_START, but also loop through each channel */
#define SIFT3D_IM_LOOP_LIMITED_START_C(im, x, y, z, c, x_start, x_end, \
			      y_start, y_end, z_start, z_end) \
        SIFT3D_IM_LOOP_LIMITED_START(im, x, y, z, x_start, x_end, \
			      y_start, y_end, z_start, z_end) \
        for ((c) = 0; (c) < (im)->nc; (c)++) {
                                

// Delimit an SIFT3D_IM_LOOP_START or SIFT3D_IM_LOOP_LIMITED_START
#define SIFT3D_IM_LOOP_END }}}

// Delimited an SIFT3D_IM_LOOP_START_C or SIFT3D_IM_LOOP_LIMITED_START_C
#define SIFT3D_IM_LOOP_END_C SIFT3D_IM_LOOP_END }

/* Check if a point is within the boundaries of an image */
#define IM_CONTAINS(im, x, y, z) \
        ((x) >= 0 && (y) >= 0 && (z) >= 0 && (x) < (im)->nx && \
         (y) < (im)->ny && (z) < (im)->nz)

/* Take the Cartesian gradient of an image at [x, y, z, c]. The voxel cannot be
 * on the boundary. */
#define SIFT3D_IM_GET_GRAD(im, x, y, z, c, vd) \
		(vd)->x = 0.5f * (SIFT3D_IM_GET_VOX(im, (x) + 1, y, z, c) - \
			   SIFT3D_IM_GET_VOX(im, (x) - 1, y, z, c)); \
		(vd)->y = 0.5f * (SIFT3D_IM_GET_VOX(im, x, (y) + 1, z, c) - \
			   SIFT3D_IM_GET_VOX(im, x, (y) - 1, z, c)); \
		(vd)->z = 0.5f * (SIFT3D_IM_GET_VOX(im, x, y, (z) + 1, c) - \
			   SIFT3D_IM_GET_VOX(im, x, y, (z) - 1, c))

/* Get the Hessian of an image at [x, y, z]. The voxel cannot be on the 
 * boundary. */
#define SIFT3D_IM_GET_HESSIAN(im, x, y, z, c, H, type) \
   /* Dxx */ \
    SIFT3D_MAT_RM_GET(H, 0, 0, type) = (type) (0.25f * \
                                (SIFT3D_IM_GET_VOX(im, x + 1, y, z, c) - \
				 2 * SIFT3D_IM_GET_VOX(im, x, y, z, c) + \
				 SIFT3D_IM_GET_VOX(im, x - 1, y, z, c))); \
    /* Dxy */ \
    SIFT3D_MAT_RM_GET(H, 0, 1, type) = (type) (0.25f * \\
                                (SIFT3D_IM_GET_VOX(im, x + 1, y + 1, z, c) - \
				 SIFT3D_IM_GET_VOX(im, x - 1, y + 1, z, c) + \
				 SIFT3D_IM_GET_VOX(im, x - 1, y - 1, z, c) - \
				 SIFT3D_IM_GET_VOX(im, x + 1, y - 1, z, c))); \
    /* Dxz */ \
    SIFT3D_MAT_RM_GET(H, 0, 2, type) = (type) (0.25f * \
                                (SIFT3D_IM_GET_VOX(im, x + 1, y, z + 1, c) - \
				 SIFT3D_IM_GET_VOX(im, x - 1, y, z + 1, c) + \
				 SIFT3D_IM_GET_VOX(im, x - 1, y, z - 1, c) - \
				 SIFT3D_IM_GET_VOX(im, x + 1, y, z - 1, c))); \
    /* Dyx */ \
    SIFT3D_MAT_RM_GET(H, 1, 0, type) = SIFT3D_MAT_RM_GET(H, 0, 1, type); \
    /* Dyy */ \
    SIFT3D_MAT_RM_GET(H, 1, 1, type) = (type) (0.25f * \
                                (SIFT3D_IM_GET_VOX(im, x, y + 1, z, c) - \
				 2 * SIFT3D_IM_GET_VOX(im, x, y, z, c) + \
				 SIFT3D_IM_GET_VOX(im, x, y - 1, z, c))); \
    /* Dyz */ \
    SIFT3D_MAT_RM_GET(H, 1, 2, type) = (type) (0.25f * \
                                (SIFT3D_IM_GET_VOX(im, x, y + 1, z + 1, c) - \
				 SIFT3D_IM_GET_VOX(im, x, y - 1, z + 1, c) + \
				 SIFT3D_IM_GET_VOX(im, x, y - 1, z - 1, c) - \
				 SIFT3D_IM_GET_VOX(im, x, y + 1, z - 1, c))); \
    /* Dzx */ \
    SIFT3D_MAT_RM_GET(H, 2, 0, type) = SIFT3D_MAT_RM_GET(H, 0, 2, type); \
    /* Dzy */ \
    SIFT3D_MAT_RM_GET(H, 2, 1, type) = SIFT3D_MAT_RM_GET(H, 1, 2, type); \
    /* Dzz */ \
    SIFT3D_MAT_RM_GET(H, 2, 2, type) = (type) (0.25f * \
                                (SIFT3D_IM_GET_VOX(im, x, y, z + 1, c) - \
				 2 * SIFT3D_IM_GET_VOX(im, x, y, z, c) + \
				 SIFT3D_IM_GET_VOX(im, x, y, z - 1, c)))

// Get a pointer to an image struct at pyramid level [o, s]
#define SIFT3D_PYR_IM_GET(pyr, o, s) ((pyr)->levels + \
						((o) - (pyr)->first_octave) * \
						(pyr)->num_levels + ((s) - (pyr)->first_level))

// Get the index of the last octave of a Pyramid struct
#define SIFT3D_PYR_LAST_OCTAVE(pyr) \
        ((pyr)->first_octave + (pyr)->num_octaves - 1)

// Get the index of the last level of a Pyramid struct
#define SIFT3D_PYR_LAST_LEVEL(pyr) \
        ((pyr)->first_level + (pyr)->num_levels - 1)

// Loop through all levels of a given pyramid
#define SIFT3D_PYR_LOOP_START(pyr, o, s) \
	for ((o) = (pyr)->first_octave; (o) <= SIFT3D_PYR_LAST_OCTAVE(pyr); \
                (o)++) { \
	for ((s) = (pyr)->first_level; (s) <= SIFT3D_PYR_LAST_LEVEL(pyr); \
                (s)++) {

// Loop from the specified (inclusive) limits of a given pyramid
#define SIFT3D_PYR_LOOP_LIMITED_START(o, s, o_start, o_end, s_start, s_end) \
	for ((o) = (o_start); (o) <= (o_end); (o)++) {	\
	for ((s) = (s_start); (s) <= (s_end); (s)++) {

// Delimit a SIFT3D_PYR_LOOP
#define SIFT3D_PYR_LOOP_END }}

// Delimit the first level of a SIFT3D_PYR_LOOP
#define SIFT3D_PYR_LOOP_SCALE_END }

// Delimit the second level of a PYR_LOOP
#define SIFT3D_PYR_LOOP_OCTAVE_END }

// Get a pointer to the incremental Gaussian filter for level s
#define SIFT3D_GAUSS_GET(gss, s) \
	((gss)->gauss_octave + (s - (gss)->first_level))

/* Resize a slab. If SIFT3D_SLAB_SIZE is defined, add
* elements in increments of that number. Otherwise,
* use a default of 500. This macro is meant to be
* used whether or not the slab buffer actually needs
* resizing -- it checks for that. */
#ifndef SIFT3D_SLAB_LEN
#define SIFT3D_SLAB_LEN 500
#endif
#define SIFT3D_RESIZE_SLAB(slab, num_new, size) \
{ \
        const size_t slabs_new = ( (size_t) (num_new) + SIFT3D_SLAB_LEN - 1) / \
                        SIFT3D_SLAB_LEN; \
        const size_t size_new = slabs_new * SIFT3D_SLAB_LEN * (size); \
\
	if (size_new != (slab)->buf_size) { \
\
		/* Re-initialize if the new size is 0 */ \
		if (size_new == 0) { \
			cleanup_Slab(slab); \
			init_Slab(slab); \
		/* Else allocate new memory */ \
		} else if (((slab)->buf = SIFT3D_safe_realloc((slab)->buf, \
			size_new)) == NULL) { \
			return SIFT3D_FAILURE; \
		} \
		(slab)->buf_size = size_new; \
	} \
	(slab)->num = (num_new); \
}

// Nested loop through all elements of a matrix
#define SIFT3D_MAT_RM_LOOP_START(mat, row, col) \
	for ((row) = 0; (row) < (mat)->num_rows; (row)++) { \
	for ((col) = 0; (col) < (mat)->num_cols; (col)++) {

// Delmit a MAT_LOOP
#define SIFT3D_MAT_RM_LOOP_END }}

// Delimit the first level of a MAT_LOOP
#define SIFT3D_MAT_RM_LOOP_COL_END }

// Delmit the second level of a MAT_LOOP
#define SIFT3D_MAT_RM_LOOP_ROW_END } 

// Get the index of an element in a dense matrix, in row-major order.
#define SIFT3D_MAT_RM_GET_IDX(mat, row, col) \
        ((col) + (row) * (mat)->num_cols)

// Get an element from a dense matrix in row-major order. Type must
// be "double", "float", or "int."
#define SIFT3D_MAT_RM_GET(mat, row, col, type) ((mat)->u.data_ ## type \
	[SIFT3D_MAT_RM_GET_IDX(mat, row, col)])

/* Execute the macro MACRO, with the first argument set to the type of mat. If
 * there is an error, goto err_label. */
#define SIFT3D_MAT_RM_TYPE_MACRO(mat, err_label, MACRO, ...) \
        switch ((mat)->type) { \
                case SIFT3D_DOUBLE: \
                        MACRO(double, ## __VA_ARGS__) \
                        break; \
                case SIFT3D_FLOAT: \
                        MACRO(float, ## __VA_ARGS__) \
                        break; \
                case SIFT3D_INT: \
                        MACRO(int, ## __VA_ARGS__) \
                        break; \
                default: \
                        SIFT3D_ERR("imutil: unknown matrix type \n"); \
                        goto err_label; \
        } \

// Convert a vector from Cartesian to Spherical coordinates.
#define SIFT3D_CVEC_TO_SVEC(cvec, svec) { \
	(svec)->mag = sqrtf((cvec)->x * (cvec)->x + (cvec)->y * (cvec)->y + \
						(cvec)->z * (cvec)->z); \
	(svec)->az = fmodf(atan2f((cvec)->y, (cvec)->x) + SIFT3D_AZ_MAX_F, \
					  SIFT3D_AZ_MAX_F); \
	(svec)->po = fmodf(acosf((cvec)->z / ((svec)->mag + FLT_EPSILON)), \
		     SIFT3D_PO_MAX_F); \
}

// Convert a vector from Spherical to Cartesian coordinates
#define SIFT3D_SVEC_TO_CVEC(svec, cvec) { \
	(cvec)->x = (svec)->mag * sinf((svec)->po) * cosf((svec)->az); \
	(cvec)->y = (svec)->mag * sinf((svec)->po) * sinf((svec)->az); \
	(cvec)->z = (svec)->mag * cosf((svec)->po); \
}

// Return the L2 norm of a Cartesian coordinate vector
#define SIFT3D_CVEC_L2_NORM(cvec) \
	sqrtf((cvec)->x * (cvec)->x + (cvec)->y * (cvec)->y + \
	(cvec)->z * (cvec)->z)

// Return the square of the  L2 norm of a Cartesian coordinate vector
#define SIFT3D_CVEC_L2_NORM_SQ(cvec) \
	((cvec)->x * (cvec)->x + (cvec)->y * (cvec)->y + \
	(cvec)->z * (cvec)->z)

// Scale a Cartesian coordinate vector by a constant factor
#define SIFT3D_CVEC_SCALE(cvec, a) { \
    (cvec)->x = (cvec)->x * a; \
    (cvec)->y = (cvec)->y * a; \
    (cvec)->z = (cvec)->z * a; \
}

// Operate element-wise on two Cartesian coordinate vectors, cc = ca op cb
#define SIFT3D_CVEC_OP(ca, cb, op, cc) { \
    (cc)->x = (ca)->x op (cb)->x; \
    (cc)->y = (ca)->y op (cb)->y; \
    (cc)->z = (ca)->z op (cb)->z; \
}

// Return the dot product of two Cartesian coordinate 
// vectors
#define SIFT3D_CVEC_DOT(in1, in2) \
	((in1)->x * (in2)->x + (in1)->y * (in2)->y + (in1)->z * (in2)->z)

// Take the cross product of two Cartesian coordinate
// vectors, as out = in1 X in2
#define SIFT3D_CVEC_CROSS(in1, in2, out) { \
	(out)->x = (in1)->y * (in2)->z - (in1)->z * (in2)->y; \
	(out)->y = (in1)->z * (in2)->x - (in1)->x * (in2)->z; \
	(out)->z = (in1)->x * (in2)->y - (in1)->y * (in2)->x; \
} 

// Evaluates to true (nonzero) if im contains cvec, false otherwise
#define SIFT3D_IM_CONTAINS_CVEC(im, cvec) ( \
        (cvec)->x >= 0 || (cvec)->y >= 0 || (cvec)->z >= 0 || \
        (cvec)->x < (float) (im)->nx || \
        (cvec)->y < (float) (im)->ny || \
        (cvec)->z < (float) (im)->nz \
)

/* Computes v_out = mat * v_in. Note that mat must be of FLOAT
 * type, since this is the only type available for vectors. 
 * Also note that mat must be (3 x 3). */
#define SIFT3D_MUL_MAT_RM_CVEC(mat, v_in, v_out) { \
	(v_out)->x = SIFT3D_MAT_RM_GET(mat, 0, 0, float) * (v_in)->x + \
	    	     SIFT3D_MAT_RM_GET(mat, 0, 1, float) * (v_in)->y + \
                     SIFT3D_MAT_RM_GET(mat, 0, 2, float) * (v_in)->z; \
	\
	(v_out)->y = SIFT3D_MAT_RM_GET(mat, 1, 0, float) * (v_in)->x + \
                     SIFT3D_MAT_RM_GET(mat, 1, 1, float) * (v_in)->y + \
                     SIFT3D_MAT_RM_GET(mat, 1, 2, float) * (v_in)->z; \
	\
	(v_out)->z = SIFT3D_MAT_RM_GET(mat, 2, 0, float) * (v_in)->x + \
                     SIFT3D_MAT_RM_GET(mat, 2, 1, float) * (v_in)->y + \
                     SIFT3D_MAT_RM_GET(mat, 2, 2, float) * (v_in)->z; \
}

#ifdef __cplusplus
}
#endif

#endif
