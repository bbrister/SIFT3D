/* macros.h
* ----------------------------------------------------------------
* Rice MRI Team
* ----------------------------------------------------------------
* This header defines preprocessor macros.
*-----------------------------------------------------------------
* Created: Blaine Rister 12/26/2013
* Last updated: Blaine Rister 11/16/2013
*/

#include "types.h"

#ifndef _MACROS_H
#define _MACROS_H

// Function return macros
#define SINGULAR 1
#define SUCCESS 0
#define FAILURE -1

// Truth macros
#define TRUE 1
#define FALSE 0

// Math macros
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define AZ_MAX_F (2 * UTIL_PI_F)
#define PO_MAX_F UTIL_PI_F

// Compiler flags
#ifdef __GNUC__
#define IGNORE_UNUSED __attribute__((unused))
#else
#warning("The internal macros of this library were not defined for " \
	 "your compiler. The code will still work, but you may see " \
	 "warnings or sub-optimal performance. \n")
#endif

// Get the index of an [x,y,z] pair in an image 
#define IM_GET_IDX(im, x, y, z) ((x) * (im)->x_stride + (y) * (im)->y_stride + \
				(z) * (im)->z_stride)

// Get the value of voxel [x,y,z] in an image 
#define IM_GET_VOX(im, x, y, z) ((im)->data[IM_GET_IDX((im), (x), (y), (z))])

// Loop through an image in x, z, y order. Delmit with IM_LOOP_END
#define IM_LOOP_START(im, x, y, z) \
	for ((z) = 0; (z) < (im)->nz; (z)++) {	\
	for ((y) = 0; (y) < (im)->ny; (y)++) {	\
	for ((x) = 0; (x) < (im)->nx; (x)++) {

/* Loop through an image iterating with the (inclusive) x, y, z bounds given.
 * Delimit with IM_LOOP_END. */
#define IM_LOOP_LIMITED_START(im, x, y, z, x_start, x_end, \
							  y_start, y_end, z_start, z_end) \
	for ((z) = z_start; (z) <= z_end; (z)++) { \
	for ((y) = y_start; (y) <= y_end; (y)++) { \
	for ((x) = x_start; (x) <= x_end; (x)++)	{		

// Delimit an IM_LOOP
#define IM_LOOP_END }}}

/* Take the Cartesian gradient of an image at [x, y, z]. The voxel cannot be on
 * the boundary. */
#define IM_GET_GRAD(im, x, y, z, vd) \
		(vd)->x = 0.5f * (IM_GET_VOX(im, x + 1, y, z) - \
			   IM_GET_VOX(im, x - 1, y, z)); \
		(vd)->y = 0.5f * (IM_GET_VOX(im, x, y + 1, z) - \
			   IM_GET_VOX(im, x, y - 1, z)); \
		(vd)->z = 0.5f * (IM_GET_VOX(im, x, y, z + 1) - \
			   IM_GET_VOX(im, x, y, z - 1))

/* Get the Hessian of an image at [x, y, z]. The voxel cannot be on the 
 * boundary. */
#define IM_GET_HESSIAN(im, x, y, z, H, type) \
   /* Dxx */ \
    MAT_RM_GET(H, 0, 0, type) = (type) (0.25f * (IM_GET_VOX(im, x + 1, y, z) - \
				 2 * IM_GET_VOX(im, x, y, z) + \
				 IM_GET_VOX(im, x - 1, y, z))); \
    /* Dxy */ \
    MAT_RM_GET(H, 0, 1, type) = (type) (0.25f * (IM_GET_VOX(im, x + 1, y + 1, z) - \
				 IM_GET_VOX(im, x - 1, y + 1, z) + \
				 IM_GET_VOX(im, x - 1, y - 1, z) - \
				 IM_GET_VOX(im, x + 1, y - 1, z))); \
    /* Dxz */ \
    MAT_RM_GET(H, 0, 2, type) = (type) (0.25f * (IM_GET_VOX(im, x + 1, y, z + 1) - \
				 IM_GET_VOX(im, x - 1, y, z + 1) + \
				 IM_GET_VOX(im, x - 1, y, z - 1) - \
				 IM_GET_VOX(im, x + 1, y, z - 1))); \
    /* Dyx */ \
    MAT_RM_GET(H, 1, 0, type) = MAT_RM_GET(H, 0, 1, type); \
    /* Dyy */ \
    MAT_RM_GET(H, 1, 1, type) = (type) (0.25f * (IM_GET_VOX(im, x, y + 1, z) - \
				 2 * IM_GET_VOX(im, x, y, z) + \
				 IM_GET_VOX(im, x, y - 1, z))); \
    /* Dyz */ \
    MAT_RM_GET(H, 1, 2, type) = (type) (0.25f * (IM_GET_VOX(im, x, y + 1, z + 1) - \
				 IM_GET_VOX(im, x, y - 1, z + 1) + \
				 IM_GET_VOX(im, x, y - 1, z - 1) - \
				 IM_GET_VOX(im, x, y + 1, z - 1))); \
    /* Dzx */ \
    MAT_RM_GET(H, 2, 0, type) = MAT_RM_GET(H, 0, 2, type); \
    /* Dzy */ \
    MAT_RM_GET(H, 2, 1, type) = MAT_RM_GET(H, 1, 2, type); \
    /* Dzz */ \
    MAT_RM_GET(H, 2, 2, type) = (type) (0.25f * (IM_GET_VOX(im, x, y, z + 1) - \
				 2 * IM_GET_VOX(im, x, y, z) + \
				 IM_GET_VOX(im, x, y, z - 1)))

// Get a pointer to an image struct at pyramid level [o, s]
#define PYR_IM_GET(pyr, o, s) ((pyr)->levels + \
						((o) - (pyr)->first_octave) * \
						(pyr)->num_levels + ((s) - (pyr)->first_level))

// Loop through all levels of a given pyramid
#define PYR_LOOP_START(pyr, o, s) \
	for ((o) = (pyr)->first_octave; (o) <= (pyr)->last_octave; (o)++) { \
	for ((s) = (pyr)->first_level; (s) <= (pyr)->last_level; (s)++) {

// Loop from the specified (inclusive) limits of a given pyramid
#define PYR_LOOP_LIMITED_START(o, s, o_start, o_end, s_start, s_end) \
	for ((o) = (o_start); (o) <= (o_end); (o)++) {	\
	for ((s) = (s_start); (s) <= (s_end); (s)++) {

// Delimit a PYR_LOOP
#define PYR_LOOP_END }}

// Delimit the first level of a PYR_LOOP
#define PYR_LOOP_SCALE_END }

// Delimit the second level of a PYR_LOOP
#define PYR_LOOP_OCTAVE_END }

// Get a pointer to the incremental Gaussian filter for level s
#define GAUSS_GET(gss, s) \
	((gss)->gauss_octave + (s - (gss)->first_level))

/* Resize a slab. If SLAB_SIZE is defined, add
* elements in increments of that number. Otherwise,
* use a default of 500. This macro is meant to be
* used whether or not the slab buffer actually needs
* resizing -- it checks for that. */
#ifndef SLAB_SIZE
#define SLAB_SIZE 500
#endif
#define RESIZE_SLAB(slab, num, size) \
	if ((num) > (slab)->buf_length) { \
		if (((slab)->buf = realloc((slab)->buf, ((slab)->buf_length + \
			SLAB_SIZE) * (size))) == NULL) \
			return FAILURE; \
		(slab)->buf_length += SLAB_SIZE; \
	} \
	(slab)->num = (num)

// As above, but can be used on a Keypoint_store
#define RESIZE_KP_STORE(store, num, size) \
	RESIZE_SLAB(&(store)->slab, num, size); \
	(store)->buf = (store)->slab.buf
#endif

// Nested loop through all elements of a matrix
#define MAT_RM_LOOP_START(mat, row, col) \
	for ((row) = 0; (row) < (mat)->num_rows; (row)++) { \
	for ((col) = 0; (col) < (mat)->num_cols; (col)++) {

// Delmit a MAT_LOOP
#define MAT_RM_LOOP_END }}

// Delimit the first level of a MAT_LOOP
#define MAT_RM_LOOP_COL_END }

// Delmit the second level of a MAT_LOOP
#define MAT_RM_LOOP_ROW_END } 

// Get an element from a dense matrix in row-major order. Type must
// be "double", "float", or "int."
#define MAT_RM_GET(mat, row, col, type) ((mat)->u.data_ ## type \
	[(col) + (row) * (mat)->num_cols])

// Convert a vector from Cartesian to Spherical coordinates.
#define CVEC_TO_SVEC(cvec, svec) { \
	(svec)->mag = sqrtf((cvec)->x * (cvec)->x + (cvec)->y * (cvec)->y + \
						(cvec)->z * (cvec)->z); \
	(svec)->az = fmodf(atan2f((cvec)->y, (cvec)->x) + AZ_MAX_F, \
					  AZ_MAX_F); \
	(svec)->po = fmodf(acosf((cvec)->z / ((svec)->mag + FLT_EPSILON)), \
		     PO_MAX_F); \
}

// Convert a vector from Spherical to Cartesian coordinates
#define SVEC_TO_CVEC(svec, cvec) { \
	(cvec)->x = (svec)->mag * sinf((svec)->po) * cosf((svec)->az); \
	(cvec)->y = (svec)->mag * sinf((svec)->po) * sinf((svec)->az); \
	(cvec)->z = (svec)->mag * cosf((svec)->po); \
}

// Return the L2 norm of a Cartesian coordinate vector
#define CVEC_L2_NORM(cvec) \
	sqrtf((cvec)->x * (cvec)->x + (cvec)->y * (cvec)->y + \
	(cvec)->z * (cvec)->z)

// Return the square of the  L2 norm of a Cartesian coordinate vector
#define CVEC_L2_NORM_SQ(cvec) \
	((cvec)->x * (cvec)->x + (cvec)->y * (cvec)->y + \
	(cvec)->z * (cvec)->z)

// Scale a Cartesian coordinate vector by a constant factor
#define CVEC_SCALE(cvec, a) \
    (cvec)->x = (cvec)->x * a; \
    (cvec)->y = (cvec)->y * a; \
    (cvec)->z = (cvec)->z * a

// Operate element-wise on two Cartesian coordinate vectors, cc = ca op cb
#define CVEC_OP(ca, cb, op, cc) \
    (cc)->x = (ca)->x op (cb)->x; \
    (cc)->y = (ca)->y op (cb)->y; \
    (cc)->z = (ca)->z op (cb)->z

// Return the dot product of two Cartesian coordinate 
// vectors
#define CVEC_DOT(in1, in2) \
	((in1)->x * (in2)->x + (in1)->y * (in2)->y + (in1)->z * (in2)->z)

// Take the cross product of two Cartesian coordinate
// vectors, as out = in1 X in2
#define CVEC_CROSS(in1, in2, out) { \
	(out)->x = (in1)->y * (in2)->z - (in1)->z * (in2)->y; \
	(out)->y = (in1)->z * (in2)->x - (in1)->x * (in2)->z; \
	(out)->z = (in1)->x * (in2)->y - (in1)->y * (in2)->x; \
} 

/* Computes v_out = mat * v_in. Note that mat must be of FLOAT
 * type, since this is the only type available for vectors. 
 * Also note that mat must be (3 x 3). */
#define MUL_MAT_RM_CVEC(mat, v_in, v_out) { \
	(v_out)->x = MAT_RM_GET(mat, 0, 0, float) * (v_in)->x + \
	    	     MAT_RM_GET(mat, 0, 1, float) * (v_in)->y + \
			     MAT_RM_GET(mat, 0, 2, float) * (v_in)->z;	\
	\
	(v_out)->y = MAT_RM_GET(mat, 1, 0, float) * (v_in)->x + \
			     MAT_RM_GET(mat, 1, 1, float) * (v_in)->y + \
			     MAT_RM_GET(mat, 1, 2, float) * (v_in)->z;	\
	\
	(v_out)->z = MAT_RM_GET(mat, 2, 0, float) * (v_in)->x + \
			     MAT_RM_GET(mat, 2, 1, float) * (v_in)->y + \
			     MAT_RM_GET(mat, 2, 2, float) * (v_in)->z; \
}

// Loop over all bins in a gradient histogram. If ICOS_HIST is defined, p
// is not referenced
#ifdef ICOS_HIST
#define HIST_LOOP_START(a, p) \
	for ((a) = 0; (a) < HIST_NUMEL; (a)++) { p = p; {
#else
#define HIST_LOOP_START(a, p) \
	for ((p) = 0; (p) < NBINS_PO; (p)++) { \
	for ((a) = 0; (a) < NBINS_AZ; (a)++) {
#endif

// Delimit a HIST_LOOP
#define HIST_LOOP_END }}

// Get an element from a gradient histogram. If ICOS_HIST is defined, p
// is not referenced
#ifdef ICOS_HIST
#define HIST_GET(hist, a, p) ((hist)->bins[a])
#else
#define HIST_GET(hist, az_bin, po_bin) ((hist)->bins[ (az_bin) + \
	(po_bin) * NBINS_AZ])
#endif
