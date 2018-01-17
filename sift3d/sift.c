/* -----------------------------------------------------------------------------
 * sift.c
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This file contains all routines needed to initialize, delete, 
 * and run the SIFT3D detector and descriptor. It also contains routines for
 * matching SIFT3D features and drawing the results.
 * -----------------------------------------------------------------------------
 */

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <getopt.h>
#include "imtypes.h"
#include "immacros.h"
#include "imutil.h"
#include "sift.h"

/* Implementation options */
//#define SIFT3D_ORI_SOLID_ANGLE_WEIGHT // Weight bins by solid angle
//#define SIFT3D_MATCH_MAX_DIST 0.3 // Maximum distance between matching features 
//#define CUBOID_EXTREMA // Search for extrema in a cuboid region

/* Internal return codes */
#define REJECT 1

/* Default SIFT3D parameters. These may be overriden by 
 * the calling appropriate functions. */
const double peak_thresh_default = 0.1; // DoG peak threshold
const int num_kp_levels_default = 3; // Number of levels per octave in which keypoints are found
const double corner_thresh_default = 0.4; // Minimum corner score
const double sigma_n_default = 1.15; // Nominal scale of input data
const double sigma0_default = 1.6; // Scale of the base octave

/* SIFT3D option names */
const char opt_peak_thresh[] = "peak_thresh";
const char opt_corner_thresh[] = "corner_thresh";
const char opt_num_kp_levels[] = "num_kp_levels";
const char opt_sigma_n[] = "sigma_n";
const char opt_sigma0[] = "sigma0";

/* Internal parameters */
const double max_eig_ratio =  0.90;	// Maximum ratio of eigenvalue magnitudes
const double ori_grad_thresh = 1E-10;   // Minimum norm of average gradient
const double bary_eps = FLT_EPSILON * 1E1;	// Error tolerance for barycentric coordinates
const double ori_sig_fctr = 1.5;        // Ratio of window parameter to keypoint scale
const double ori_rad_fctr =  3.0; // Ratio of window radius to parameter
const double desc_sig_fctr = 7.071067812; // See ori_sig_fctr, 5 * sqrt(2)
const double desc_rad_fctr = 2.0;  // See ori_rad_fctr
const double trunc_thresh = 0.2f * 128.0f / DESC_NUMEL; // Descriptor truncation threshold

/* Internal math constants */
const double gr = 1.6180339887; // Golden ratio

/* Get the index of bin j from triangle i */
#define MESH_GET_IDX(mesh, i, j) \
	((mesh)->tri[i].idx[j])

/* Get bin j from triangle i */
#define MESH_HIST_GET(mesh, hist, i, j) \
	((hist)->bins[MESH_GET_IDX(mesh, i, j)])

/* Clamp out of bounds polar accesses to the first or last element.
 * Note that the polar histogram is NOT circular. */
#define HIST_GET_PO(hist, a, p) \
			 ((p) < 0 ? \
			 HIST_GET(hist, ((a) + NBINS_AZ / 2) % NBINS_AZ, 1) : \
			 (p) >= NBINS_PO ? \
		 	 HIST_GET(hist, ((a) + NBINS_AZ / 2) % NBINS_AZ, \
			 NBINS_PO - 1) : \
		     HIST_GET(hist, a, p))

/* Convert out of bounds azimuthal accesses circularly, e.g. -1 goes
 * to NBINS_AZ - 1, NBINS_AZ goes to 0. This algorithm does not work
 * if the indices wrap around more than once. */
#define HIST_GET_AZ(hist, a, p)	\
			 HIST_GET_PO(hist, ((a) + NBINS_AZ) % NBINS_AZ, p)

/* Loop through a spherical image region. im and [x, y, z] are defined as
 * above. vcenter is a pointer to a Cvec specifying the center of the window.
 * rad is the radius of the window. vdisp is a pointer to a Cvec storing
 * the displacement from the window center. sqdisp is a float storing the
 * squared Euclidean distance from the window center.
 *
 * Note that the sphere is defined in real-world coordinates, i.e. those
 * with units (1, 1, 1). Thus, rad, sq_dist, and vdisp are defined in these
 * coordinates as well. However, x, y, z, and vcenter are defined in image
 * space.
 *
 * Delimit with IM_LOOP_SPHERE_END. */
#define IM_LOOP_SPHERE_START(im, x, y, z, vcenter, rad, vdisp, sq_dist) \
{ \
        const float uxf = (float) (im)->ux; \
        const float uyf = (float) (im)->uy; \
        const float uzf = (float) (im)->uz; \
	const int x_start = SIFT3D_MAX(floorf((vcenter)->x - (rad) / uxf), 1); \
	const int x_end   = SIFT3D_MIN(ceilf((vcenter)->x + (rad) / uxf),  \
                im->nx - 2); \
	const int y_start = SIFT3D_MAX(floorf((vcenter)->y - (rad) / uyf), 1); \
	const int y_end   = SIFT3D_MIN(ceilf((vcenter)->y + (rad) / uyf), \
                im->ny - 2); \
	const int z_start = SIFT3D_MAX(floorf((vcenter)->z - (rad) / uzf), 1); \
	const int z_end   = SIFT3D_MIN(ceilf((vcenter)->z + (rad) / uzf), \
                im->nz - 2); \
	SIFT3D_IM_LOOP_LIMITED_START(im, x, y, z, x_start, x_end, y_start, \
                 y_end, z_start, z_end) \
                (vdisp)->x = ((float) x - (vcenter)->x) * uxf; \
                (vdisp)->y = ((float) y - (vcenter)->y) * uyf; \
                (vdisp)->z = ((float) z - (vcenter)->z) * uzf; \
                (sq_dist) = SIFT3D_CVEC_L2_NORM_SQ(vdisp); \
                if ((sq_dist) > (rad) * (rad)) \
	                continue; \

#define IM_LOOP_SPHERE_END SIFT3D_IM_LOOP_END }

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
#define HIST_GET_IDX(a, p) (a)
#else
#define HIST_GET_IDX(a, p) ((a) + (p) * NBINS_AZ)
#endif
#define HIST_GET(hist, a, p) ((hist)->bins[HIST_GET_IDX(a, p)])

// Get a column index in the matrix representation of a 
// SIFT3D_Descriptor_store struct
#define DESC_MAT_GET_COL(hist_idx, a, p) \
        (((hist_idx) * HIST_NUMEL) + HIST_GET_IDX(a, p) + IM_NDIMS)

// As SIFT3D_IM_GET_GRAD, but with physical units (1, 1, 1)
#define IM_GET_GRAD_ISO(im, x, y, z, c, vd) { \
        SIFT3D_IM_GET_GRAD(im, x, y, z, c, vd); \
        (vd)->x *=  1.0f / (float) (im)->ux; \
        (vd)->y *= 1.0f / (float) (im)->uy; \
        (vd)->z *= 1.0f / (float) (im)->uz; \
}

/* Global variables */
extern CL_data cl_data;

/* Helper routines */
static int init_geometry(SIFT3D *sift3d);
static int set_im_SIFT3D(SIFT3D *const sift3d, const Image *const im);
static int set_scales_SIFT3D(SIFT3D *const sift3d, const double sigma0,
        const double sigma_n);
static int resize_SIFT3D(SIFT3D *const sift3d, const int num_kp_levels);
static int build_gpyr(SIFT3D *sift3d);
static int build_dog(SIFT3D *dog);
static int detect_extrema(SIFT3D *sift3d, Keypoint_store *kp);
static int assign_orientations(SIFT3D *const sift3d, Keypoint_store *const kp);
static int assign_orientation_thresh(const Image *const im, 
        const Cvec *const vcenter, const double sigma, const double thresh,
        Mat_rm *const R);
static int assign_eig_ori(const Image *const im, const Cvec *const vcenter,
                          const double sigma, Mat_rm *const R, 
                          double *const conf);
static int Cvec_to_sbins(const Cvec * const vd, Svec * const bins);
static void refine_Hist(Hist *hist);
static int init_cl_SIFT3D(SIFT3D *sift3d);
static int cart2bary(const Cvec * const cart, const Tri * const tri, 
		      Cvec * const bary, float * const k);
static int scale_Keypoint(const Keypoint *const src, 
        const double *const factors, Keypoint *const dst);
static int smooth_scale_raw_input(const SIFT3D *const sift3d, const Image *const src,
        Image *const dst);
static int verify_keys(const Keypoint_store *const kp, const Image *const im);
static int keypoint2base(const Keypoint *const src, Keypoint *const dst);
static int _SIFT3D_extract_descriptors(SIFT3D *const sift3d, 
        const Pyramid *const gpyr, const Keypoint_store *const kp, 
        SIFT3D_Descriptor_store *const desc);
static void SIFT3D_desc_acc_interp(const SIFT3D * const sift3d, 
				   const Cvec * const vbins, 
				   const Cvec * const grad,
				   SIFT3D_Descriptor * const desc);
static int extract_descrip(SIFT3D *const sift3d, const Image *const im,
	   const Keypoint *const key, SIFT3D_Descriptor *const desc);
static int argv_remove(const int argc, char **argv, 
                        const unsigned char *processed);
static int extract_dense_descriptors_no_rotate(SIFT3D *const sift3d,
        const Image *const in, Image *const desc);
static int extract_dense_descriptors_rotate(SIFT3D *const sift3d,
        const Image *const in, Image *const desc);
static int extract_dense_descrip_rotate(SIFT3D *const sift3d, 
           const Image *const im, const Cvec *const vcenter, 
           const double sigma, const Mat_rm *const R, Hist *const hist);
static void vox2hist(const Image *const im, const int x, const int y,
        const int z, Hist *const hist);
static void hist2vox(Hist *const hist, const Image *const im, const int x, 
        const int y, const int z);
static int match_desc(const SIFT3D_Descriptor *const desc,
        const SIFT3D_Descriptor_store *const store, const float nn_thresh);
static int resize_SIFT3D_Descriptor_store(SIFT3D_Descriptor_store *const desc,
        const int num);

/* Initialize geometry tables. */
static int init_geometry(SIFT3D *sift3d) {

	Mat_rm V, F;
	Cvec temp1, temp2, temp3, n;
	float mag;
	int i, j;

	Mesh * const mesh = &sift3d->mesh;

	/* Verices of a regular icosahedron inscribed in the unit sphere. */
	const float vert[] = {  0,  1,  gr,
			        0, -1,  gr,
			        0,  1, -gr,
			        0, -1, -gr,
			        1,  gr,  0,
			       -1,  gr,  0,
			        1, -gr,  0,
			       -1, -gr,  0,
			       gr,   0,  1,
			      -gr,   0,  1,
			       gr,   0, -1, 
			      -gr,   0, -1 }; 

	/* Vertex triplets forming the faces of the icosahedron. */
	const float faces[] = {0, 1, 8,
    			       0, 8, 4,
    			       0, 4, 5,
    			       0, 5, 9,
    			       0, 9, 1,
    			       1, 6, 8,
			       8, 6, 10,
			       8, 10, 4,
			       4, 10, 2,
			       4, 2, 5,
			       5, 2, 11,
			       5, 11, 9,
			       9, 11, 7,
			       9, 7, 1,
			       1, 7, 6,
			       3, 6, 7,
			       3, 7, 11,
			       3, 11, 2,
			       3, 2, 10,
			       3, 10, 6};

	// Initialize matrices
	if (init_Mat_rm_p(&V, vert, ICOS_NVERT, 3, SIFT3D_FLOAT, 
		SIFT3D_FALSE) ||
	    init_Mat_rm_p(&F, faces, ICOS_NFACES, 3, SIFT3D_FLOAT, 
	    	SIFT3D_FALSE))
		return SIFT3D_FAILURE;
			    
	// Initialize triangle memory
        init_Mesh(mesh);
	if ((mesh->tri = (Tri *) SIFT3D_safe_realloc(mesh->tri, 
                ICOS_NFACES * sizeof(Tri))) == NULL)
		return SIFT3D_FAILURE;
 
	// Populate the triangle struct for each face
	for (i = 0; i < ICOS_NFACES; i++) {

		Tri * const tri = mesh->tri + i;	
		Cvec * const v = tri->v;

		// Initialize the vertices
		for (j = 0; j < 3; j++) {

			const float mag_expected = sqrt(1 + gr * gr);
			int * const idx = tri->idx + j;

			*idx = SIFT3D_MAT_RM_GET(&F, i, j, float);

			// Initialize the vector
			v[j].x = SIFT3D_MAT_RM_GET(&V, *idx, 0, float);
			v[j].y = SIFT3D_MAT_RM_GET(&V, *idx, 1, float);
			v[j].z = SIFT3D_MAT_RM_GET(&V, *idx, 2, float);

			// Normalize to unit length
			mag = SIFT3D_CVEC_L2_NORM(v + j);
			assert(fabsf(mag - mag_expected) < 1E-10);
			SIFT3D_CVEC_SCALE(v + j, 1.0f / mag);
		}

		// Compute the normal vector at v[0] as  (V2 - V1) X (V1 - V0)
		SIFT3D_CVEC_OP(v + 2, v + 1, -, &temp1);
		SIFT3D_CVEC_OP(v + 1, v, -, &temp2);
		SIFT3D_CVEC_CROSS(&temp1, &temp2, &n);

		// Ensure this vector is facing outward from the origin
		if (SIFT3D_CVEC_DOT(&n, v) < 0) {
			// Swap two vertices
			temp1 = v[0];
			v[0] = v[1];
			v[1] = temp1;

			// Compute the normal again
			SIFT3D_CVEC_OP(v + 2, v + 1, -, &temp1);
			SIFT3D_CVEC_OP(v + 1, v, -, &temp2);
			SIFT3D_CVEC_CROSS(&temp1, &temp2, &n);
		}
		assert(SIFT3D_CVEC_DOT(&n, v) >= 0);

		// Ensure the triangle is equilateral
		SIFT3D_CVEC_OP(v + 2, v, -, &temp3);
		assert(fabsf(SIFT3D_CVEC_L2_NORM(&temp1) - 
                        SIFT3D_CVEC_L2_NORM(&temp2)) < 1E-10);
		assert(fabsf(SIFT3D_CVEC_L2_NORM(&temp1) - 
                        SIFT3D_CVEC_L2_NORM(&temp3)) < 1E-10);
	}	
	
	return SIFT3D_SUCCESS;
}

/* Convert Cartesian coordinates to barycentric. bary is set to all zeros if
 * the problem is unstable. 
 *
 * The output value k is the constant by which the ray is multiplied to
 * intersect the supporting plane of the triangle.
 *
 * This code uses the Moller-Trumbore algorithm. */
static int cart2bary(const Cvec * const cart, const Tri * const tri, 
		      Cvec * const bary, float * const k) {

	Cvec e1, e2, t, p, q;
	float det, det_inv;

	const Cvec * const v = tri->v;

	SIFT3D_CVEC_OP(v + 1, v, -, &e1);
	SIFT3D_CVEC_OP(v + 2, v, -, &e2);
	SIFT3D_CVEC_CROSS(cart, &e2, &p);
	det = SIFT3D_CVEC_DOT(&e1, &p);

	// Reject unstable points
	if (fabsf(det) < bary_eps) {
		return SIFT3D_FAILURE;
	}

	det_inv = 1.0f / det;

	t = v[0];
	SIFT3D_CVEC_SCALE(&t, -1.0f);	

	SIFT3D_CVEC_CROSS(&t, &e1, &q);

	bary->y = det_inv * SIFT3D_CVEC_DOT(&t, &p);	
	bary->z = det_inv * SIFT3D_CVEC_DOT(cart, &q);
	bary->x = 1.0f - bary->y - bary->z;

	*k = SIFT3D_CVEC_DOT(&e2, &q) * det_inv;

#ifndef NDEBUG
	Cvec temp1, temp2, temp3;
        double residual;

        if (isnan(bary->x) || isnan(bary->y) || isnan(bary->z)) {
                printf("cart2bary: invalid bary (%f, %f, %f)\n", bary->x, 
                        bary->y, bary->z);
                //exit(1);
        }

	// Verify k * c = bary->x * v1 + bary->y * v2 + bary->z * v3
	temp1 = v[0];
	temp2 = v[1];
	temp3 = v[2];
	SIFT3D_CVEC_SCALE(&temp1, bary->x);
	SIFT3D_CVEC_SCALE(&temp2, bary->y);	
	SIFT3D_CVEC_SCALE(&temp3, bary->z);	
	SIFT3D_CVEC_OP(&temp1, &temp2, +, &temp1);
	SIFT3D_CVEC_OP(&temp1, &temp3, +, &temp1);
	SIFT3D_CVEC_SCALE(&temp1, 1.0f / *k);
	SIFT3D_CVEC_OP(&temp1, cart, -, &temp1);
        residual = SIFT3D_CVEC_L2_NORM(&temp1);
	if (residual > bary_eps) {
                printf("cart2bary: residual: %f\n", residual);
                exit(1);
        }
#endif
	return SIFT3D_SUCCESS;
}

/* Initialize a Keypoint_store for first use.
 * This does not need to be called to reuse the store
 * for a new image. */
void init_Keypoint_store(Keypoint_store *const kp) {
	init_Slab(&kp->slab);
	kp->buf = (Keypoint *) kp->slab.buf;
}

/* Initialize a Keypoint struct for use. This sets up the internal pointers,
 * and nothing else. If called on a valid Keypoint struct, it has no effect. */
int init_Keypoint(Keypoint *const key) {
        // Initialize the orientation matrix with static memory
        return init_Mat_rm_p(&key->R, key->r_data, IM_NDIMS, IM_NDIMS, 
		SIFT3D_FLOAT, SIFT3D_FALSE);
}

/* Make room for at least num Keypoint structs in kp. 
 * 
 * Note: This function must re-initialize some internal data if it was moved. 
 * This does not affect the end user, but it affects the implementation of 
 * init_Keypoint. */
int resize_Keypoint_store(Keypoint_store *const kp, const size_t num) {

        void *const buf_old = kp->slab.buf;

        // Resize the internal memory
	SIFT3D_RESIZE_SLAB(&kp->slab, num, sizeof(struct _Keypoint));
	kp->buf = kp->slab.buf; 

        // If the size has changed, re-initialize the keypoints
        if (buf_old != kp->slab.buf) { 
                int i; 
                for (i = 0; i < kp->slab.num; i++) { 
                        Keypoint *const key = kp->buf + i; 
                        if (init_Keypoint(key)) 
                                return SIFT3D_FAILURE; 
                } 
        } 

        return SIFT3D_SUCCESS;
}

/* Copy one Keypoint struct into another. */
int copy_Keypoint(const Keypoint *const src, Keypoint *const dst) {

        // Copy the shallow data 
        dst->xd = src->xd;
        dst->yd = src->yd;
        dst->zd = src->zd;
        dst->sd = src->sd;
        dst->o = src->o;
        dst->s = src->s;

        // Copy the orienation matrix
        return copy_Mat_rm(&src->R, &dst->R);
}

/* Free all memory associated with a Keypoint_store. kp cannot be
 * used after calling this function, unless re-initialized. */
void cleanup_Keypoint_store(Keypoint_store *const kp) {
        cleanup_Slab(&kp->slab);
}

/* Initialize a SIFT_Descriptor_store for first use.
 * This does not need to be called to reuse the store
 * for a new image. */
void init_SIFT3D_Descriptor_store(SIFT3D_Descriptor_store *const desc) {
	desc->buf = NULL;
}

/* Free all memory associated with a SIFT3D_Descriptor_store. desc
 * cannot be used after calling this function, unless re-initialized. */
void cleanup_SIFT3D_Descriptor_store(SIFT3D_Descriptor_store *const desc) {
        free(desc->buf);
}

/* Resize a SIFT3D_Descriptor_store to hold n descriptors. Must be initialized
 * prior to calling this function. num must be positive.
 *
 * Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
static int resize_SIFT3D_Descriptor_store(SIFT3D_Descriptor_store *const desc,
        const int num) {

        if (num < 1) {
                SIFT3D_ERR("resize_SIFT3D_Descriptor_store: invalid size: %d",
                        num);
                return SIFT3D_FAILURE;
        }

	if ((desc->buf = (SIFT3D_Descriptor *) SIFT3D_safe_realloc(desc->buf, 
		num * sizeof(SIFT3D_Descriptor))) == NULL)
                return SIFT3D_FAILURE;

	desc->num = num;
        return SIFT3D_SUCCESS;
}

/* Initializes the OpenCL data for this SIFT3D struct. This
 * increments the reference counts for shared data. */
static int init_cl_SIFT3D(SIFT3D *sift3d) {
#ifdef SIFT3D_USE_OPENCL
	cl_image_format image_format;

	// Initialize basic OpenCL platform and context info
	image_format.image_channel_order = CL_R;
	image_format.image_channel_data_type = CL_FLOAT;
	if (init_cl(&cl_data, PLATFORM_NAME_NVIDIA, CL_DEVICE_TYPE_GPU,
 		    CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, 
                    image_format))
		return SIFT3D_FAILURE;

	// Load and compile the downsampling kernel

#endif
	return SIFT3D_SUCCESS;
}

/* Sets the peak threshold, checking that it is in the interval (0, inf) */
int set_peak_thresh_SIFT3D(SIFT3D *const sift3d,
                                const double peak_thresh) {
        if (peak_thresh <= 0.0 || peak_thresh > 1) {
                SIFT3D_ERR("SIFT3D peak_thresh must be in the interval (0, 1]. "
                        "Provided: %f \n", peak_thresh);
                return SIFT3D_FAILURE;
        }

        sift3d->peak_thresh = peak_thresh;
        return SIFT3D_SUCCESS;
}

/* Sets the corner threshold, checking that it is in the interval [0, 1]. */
int set_corner_thresh_SIFT3D(SIFT3D *const sift3d,
                                const double corner_thresh) {

        if (corner_thresh < 0.0 || corner_thresh > 1.0) {
                SIFT3D_ERR("SIFT3D corner_thresh must be in the interval "
                        "[0, 1]. Provided: %f \n", corner_thresh);
                return SIFT3D_FAILURE;
        }

        sift3d->corner_thresh = corner_thresh;
        return SIFT3D_SUCCESS;
}

/* Sets the number of levels per octave. This function will resize the
 * internal data. */
int set_num_kp_levels_SIFT3D(SIFT3D *const sift3d,
                                const unsigned int num_kp_levels) {

        const Pyramid *const gpyr = &sift3d->gpyr;

        return resize_SIFT3D(sift3d, num_kp_levels);
}

/* Sets the nominal scale parameter of the input data, checking that it is 
 * nonnegative. */
int set_sigma_n_SIFT3D(SIFT3D *const sift3d,
                                const double sigma_n) {

        const double sigma0 = sift3d->gpyr.sigma0;

        if (sigma_n < 0.0) {
                SIFT3D_ERR("SIFT3D sigma_n must be nonnegative. Provided: "
                        "%f \n", sigma_n);
                return SIFT3D_FAILURE;
        }

        return set_scales_SIFT3D(sift3d, sigma0, sigma_n);
}

/* Sets the scale parameter of the first level of octave 0, checking that it
 * is nonnegative. */
int set_sigma0_SIFT3D(SIFT3D *const sift3d,
                                const double sigma0) {

        const double sigma_n = sift3d->gpyr.sigma_n;

        if (sigma0 < 0.0) {
                SIFT3D_ERR("SIFT3D sigma0 must be nonnegative. Provided: "
                        "%f \n", sigma0);
                return SIFT3D_FAILURE; 
        } 

        return set_scales_SIFT3D(sift3d, sigma0, sigma_n);
}

/* Initialize a SIFT3D struct with the default parameters. */
int init_SIFT3D(SIFT3D *sift3d) {

        Pyramid *const dog = &sift3d->dog;
        Pyramid *const gpyr = &sift3d->gpyr;
        GSS_filters *const gss = &sift3d->gss;

	// Initialize to defaults
	const double peak_thresh = peak_thresh_default;
	const double corner_thresh = corner_thresh_default;
	const int num_kp_levels = num_kp_levels_default;
	const double sigma_n = sigma_n_default;
	const double sigma0 = sigma0_default;
        const int dense_rotate = SIFT3D_FALSE;

	// First-time pyramid initialization
        init_Pyramid(dog);
        init_Pyramid(gpyr);

        // First-time filter initialization
        init_GSS_filters(gss);

        // Intialize the geometry tables
	if (init_geometry(sift3d))
		return SIFT3D_FAILURE;

	// init static OpenCL programs and contexts, if support is enabled
	if (init_cl_SIFT3D(sift3d))
		return SIFT3D_FAILURE;

	// Initialize the image data
	init_im(&sift3d->im);

	// Save data
	dog->first_level = gpyr->first_level = -1;
        sift3d->dense_rotate = dense_rotate;
        if (set_sigma_n_SIFT3D(sift3d, sigma_n) ||
                set_sigma0_SIFT3D(sift3d, sigma0) ||
                set_peak_thresh_SIFT3D(sift3d, peak_thresh) ||
                set_corner_thresh_SIFT3D(sift3d, corner_thresh) ||
                set_num_kp_levels_SIFT3D(sift3d, num_kp_levels))
                return SIFT3D_FAILURE;

	return SIFT3D_SUCCESS;
}

/* Make a deep copy of a SIFT3D struct, including all internal images. */
int copy_SIFT3D(const SIFT3D *const src, SIFT3D *const dst) {

        // Free and re-initialize dst
        cleanup_SIFT3D(dst);
        if (init_SIFT3D(dst))
                return SIFT3D_FAILURE;

        // Copy the parameters
        set_sigma_n_SIFT3D(dst, src->gpyr.sigma_n); 
        set_sigma0_SIFT3D(dst, src->gpyr.sigma0);
        if (set_peak_thresh_SIFT3D(dst, src->peak_thresh) ||
            set_corner_thresh_SIFT3D(dst, src->corner_thresh) ||
            set_num_kp_levels_SIFT3D(dst, src->gpyr.num_kp_levels))
                return SIFT3D_FAILURE;
        dst->dense_rotate = src->dense_rotate;

        // Copy the image, if any
        if (src->im.data != NULL && set_im_SIFT3D(dst, &src->im))
                return SIFT3D_FAILURE;

        // Copy the pyramids, if any
        if (copy_Pyramid(&src->gpyr, &dst->gpyr) ||
            copy_Pyramid(&src->dog, &dst->dog))
                return SIFT3D_FAILURE;

        return SIFT3D_SUCCESS;
}

/* Free all memory associated with a SIFT3D struct. sift3d cannot be reused
 * unless it is reinitialized. */
void cleanup_SIFT3D(SIFT3D *const sift3d) {

	// Clean up the image copy
	im_free(&sift3d->im);

        // Clean up the pyramids
        cleanup_Pyramid(&sift3d->gpyr);
        cleanup_Pyramid(&sift3d->dog);

        // Clean up the GSS filters
        cleanup_GSS_filters(&sift3d->gss);

        // Clean up the triangle mesh 
        cleanup_Mesh(&sift3d->mesh);

#ifdef USE_OPENCL
        // Clean up the OpenCL kernels
        cleanup_SIFT3D_cl_kernels(&sift3d->kernels);
#endif
}

/* Helper function to remove the processed arguments from argv. 
 * Returns the number of remaining arguments. */
static int argv_remove(const int argc, char **argv, 
                        const unsigned char *processed) {

        int i, new_pos;

        // Remove the processed arguments in-place
        new_pos = 0;
        for (i = 0; i < argc; i++) {
                // Skip processed arguments
                if (processed[i])
                        continue;

                // Add the unprocessed arguments to the new argv
                argv[new_pos++] = argv[i];                 
        }

        return new_pos;
}

/* Print the options for a SIFT3D struct to stdout. */
void print_opts_SIFT3D(void) {

        printf("SIFT3D Options: \n"
               " --%s [value] \n"
               "    The smallest allowed absolute DoG value, as a fraction \n"
               "        of the largest. Must be on the interval (0, 1]. \n"
               "        (default: %.2f) \n" 
               " --%s [value] \n"
               "    The smallest allowed corner score, on the interval \n"
               "        [0, 1]. (default: %.2f) \n"
               " --%s [value] \n"
               "    The number of pyramid levels per octave in which \n"
               "        keypoints are found. Must be a positive integer. \n"
               "        (default: %d) \n"
               " --%s [value] \n"
               "    The nominal scale parameter of the input data, on the \n"
               "        interval (0, inf). (default: %.2f) \n"
               " --%s [value] \n"
               "    The scale parameter of the first level of octave 0, on \n"
               "        the interval (0, inf). (default: %.2f) \n",
               opt_peak_thresh, peak_thresh_default,
               opt_corner_thresh, corner_thresh_default,
               opt_num_kp_levels, num_kp_levels_default,
               opt_sigma_n, sigma_n_default,
               opt_sigma0, sigma0_default);

}

/* Set the parameters of a SIFT3D struct from the given command line 
 * arguments. The argument SIFT3D must be initialized with
 * init_SIFT3D prior to calling this function. 
 *
 * On return, all processed SIFT3D options will be removed from argv.
 * Use argc_ret to get the number of remaining options.
 *
 * Options:
 * --peak_thresh	 - threshold on DoG extrema magnitude (double)
 * --corner_thresh - threshold on edge score (double)
 * --num_kp_levels - number of levels per octave for keypoint
 *    			candidates (int)
 * --sigma_n - base level of blurring assumed in data (double)
 * --sigma0 - level to blur base of pyramid (double)
 *
 * Parameters:
 *      argc - The number of arguments
 *      argv - An array of strings of arguments. All unproccesed arguments are
 *              permuted to the end.
 *      sift3d - The struct to be initialized
 *      check_err - If nonzero, report unrecognized options as errors
 *
 * Return value: 
 *       Returns the new number of arguments in argv, or -1 on failure. */
int parse_args_SIFT3D(SIFT3D *const sift3d,
        const int argc, char **argv, const int check_err) {

        unsigned char *processed;
        double dval;
        int c, err, ival, argc_new;

#define PEAK_THRESH 'a'
#define CORNER_THRESH 'b'
#define NUM_KP_LEVELS 'c'
#define SIGMA_N 'd'
#define SIGMA0 'e'

        // Options
        const struct option longopts[] = {
                {opt_peak_thresh, required_argument, NULL, PEAK_THRESH},
                {opt_corner_thresh, required_argument, NULL, CORNER_THRESH},
                {opt_num_kp_levels, required_argument, NULL, NUM_KP_LEVELS},
                {opt_sigma_n, required_argument, NULL, SIGMA_N},
                {opt_sigma0, required_argument, NULL, SIGMA0},
                {0, 0, 0, 0}
        };

        // Starting getopt variables 
        const int opterr_start = opterr;

        // Set the error checking behavior
        opterr = check_err;

        // Intialize intermediate data
        if ((processed = calloc(argc, sizeof(char *))) == NULL) {
                SIFT3D_ERR("parse_args_SIFT3D: out of memory \n");
                return -1;
        }
        err = SIFT3D_FALSE;

        // Process the arguments
        while ((c = getopt_long(argc, argv, "-", longopts, NULL)) != -1) {

                const int idx = optind - 1;

                // Convert the value to double and integer
                if (optarg != NULL) {
                        dval = atof(optarg);
                        ival = atoi(optarg);
                }

                switch (c) {
                        case PEAK_THRESH:
                                if (set_peak_thresh_SIFT3D(sift3d, 
                                        dval))
                                        goto parse_args_quit;

                                processed[idx - 1] = SIFT3D_TRUE;
                                processed[idx] = SIFT3D_TRUE;
                                break;
                        case CORNER_THRESH:
                                if (set_corner_thresh_SIFT3D(sift3d, 
                                        dval))
                                        goto parse_args_quit;

                                processed[idx - 1] = SIFT3D_TRUE;
                                processed[idx] = SIFT3D_TRUE;
                                break;
                        case NUM_KP_LEVELS:
                                // Check for errors                        
                                if (ival <= 0) {
                                        SIFT3D_ERR("SIFT3D num_kp_levels "
                                                "must be positive. Provided: "
                                                "%d \n", ival);
                                        goto parse_args_quit;
                                }

                                if (set_num_kp_levels_SIFT3D(sift3d, 
                                        ival))
                                        goto parse_args_quit;

                                processed[idx - 1] = SIFT3D_TRUE;
                                processed[idx] = SIFT3D_TRUE;
                                break;
                        case SIGMA_N:
                                set_sigma_n_SIFT3D(sift3d, dval);
                                processed[idx - 1] = SIFT3D_TRUE;
                                processed[idx] = SIFT3D_TRUE;
                                break;
                        case SIGMA0:
                                set_sigma0_SIFT3D(sift3d, dval);
                                processed[idx - 1] = SIFT3D_TRUE;
                                processed[idx] = SIFT3D_TRUE;
                                break;
                        case '?':
                        default:
                                if (!check_err)
                                        continue;
                                err = SIFT3D_TRUE;
                }
        }

#undef PEAK_THRESH
#undef CORNER_THRESH
#undef NUM_KP_LEVELS
#undef SIGMA_N
#undef SIGMA0

        // Put all unprocessed options at the end
        argc_new = argv_remove(argc, argv, processed);

        // Return to the default settings
        opterr = opterr_start;

        // Clean up
        free(processed);

        // Return an error, if error checking is enabled
        if (check_err && err)
                return -1;
        
        // Reset the state of getopt
        optind = 0;

        return argc_new;

parse_args_quit:
        free(processed);
        return -1;
}

/* Helper routine to begin processing a new image. If the dimensions differ
 * from the last one, this function resizes the SIFT3D struct. */
static int set_im_SIFT3D(SIFT3D *const sift3d, const Image *const im) {

        int dims_old[IM_NDIMS];
        int i;

	const float *const data_old = sift3d->im.data;
        const Pyramid *const gpyr = &sift3d->gpyr;
        const int first_octave = sift3d->gpyr.first_octave;
        const int num_kp_levels = gpyr->num_kp_levels;

        // Make a temporary copy the previous image dimensions
        for (i = 0; i < IM_NDIMS; i++) {
                dims_old[i] = SIFT3D_IM_GET_DIMS(&sift3d->im)[i];
        }

        // Make a copy of the input image
        if (im_copy_data(im, &sift3d->im))
                return SIFT3D_FAILURE;

        // Scale the input image to [-1, 1]
        im_scale(&sift3d->im);

        // Resize the internal data, if necessary
        if ((data_old == NULL || 
                memcmp(dims_old, SIFT3D_IM_GET_DIMS(&sift3d->im), 
                        IM_NDIMS * sizeof(int))) &&
                resize_SIFT3D(sift3d, num_kp_levels))
                return SIFT3D_FAILURE;

        return SIFT3D_SUCCESS;
}

/* Helper function to set the scale parameters for a SIFT3D struct. */
static int set_scales_SIFT3D(SIFT3D *const sift3d, const double sigma0,
        const double sigma_n) {

        Pyramid *const gpyr = &sift3d->gpyr;
        Pyramid *const dog = &sift3d->dog;
        GSS_filters *const gss = &sift3d->gss;

        // Set the scales for the GSS and DOG pyramids
        if (set_scales_Pyramid(sigma0, sigma_n, gpyr) ||
                set_scales_Pyramid(sigma0, sigma_n, dog))
                return SIFT3D_FAILURE;

        // Do nothing more if we have no image
        if (sift3d->im.data == NULL)
                return SIFT3D_SUCCESS;

        // Recompute the filters
	return make_gss(gss, gpyr);
}

/* Resize a SIFT3D struct, allocating temporary storage and recompiling the 
 * filters. Does nothing unless set_im_SIFT3D was previously called. */
static int resize_SIFT3D(SIFT3D *const sift3d, const int num_kp_levels) {

        int num_octaves; 

        const Image *const im = &sift3d->im;
        Pyramid *const gpyr = &sift3d->gpyr;
        Pyramid *const dog = &sift3d->dog;
	const unsigned int num_dog_levels = num_kp_levels + 2;
	const unsigned int num_gpyr_levels = num_dog_levels + 1;
        const int first_octave = 0;
        const int first_level = -1;

	// Compute the meximum allowed number of octaves
	if (im->data != NULL) {
                // The minimum size of a pyramid level is 8 in any dimension
		const int last_octave = 
                        (int) log2((double) SIFT3D_MIN(SIFT3D_MIN(im->nx, im->ny), 
                        im->nz)) - 3 - first_octave;

                // Verify octave parameters
                if (last_octave < first_octave) {
                        SIFT3D_ERR("resize_SIFT3D: input image is too small: "
                                "must have at least 8 voxels in each "
                                "dimension \n");
                        return SIFT3D_FAILURE;
                }

                num_octaves = last_octave - first_octave + 1;
	} else {
                num_octaves = 0;
        }

	// Resize the pyramid
	if (resize_Pyramid(im, first_level, num_kp_levels,
                num_gpyr_levels, first_octave, num_octaves, gpyr) ||
	        resize_Pyramid(im, first_level, num_kp_levels, 
                num_dog_levels, first_octave, num_octaves, dog))
		return SIFT3D_FAILURE;

        // Do nothing more if we have no image
        if (im->data == NULL)
                return SIFT3D_SUCCESS;

	// Compute the Gaussian filters
	if (make_gss(&sift3d->gss, &sift3d->gpyr))
		return SIFT3D_FAILURE;

	return SIFT3D_SUCCESS;
}

/* Build the GSS pyramid on a single CPU thread */
static int build_gpyr(SIFT3D *sift3d) {

        const Image *prev;
	Sep_FIR_filter *f;
	Image *cur;
	int o, s;

	Pyramid *const gpyr = &sift3d->gpyr;
	const GSS_filters *const gss = &sift3d->gss;
	const int s_start = gpyr->first_level + 1;
	const int s_end = SIFT3D_PYR_LAST_LEVEL(gpyr);
	const int o_start = gpyr->first_octave;
	const int o_end = SIFT3D_PYR_LAST_OCTAVE(gpyr);
        const double unit = 1.0;

	// Build the first image
	cur = SIFT3D_PYR_IM_GET(gpyr, o_start, s_start - 1);
	prev = &sift3d->im;
#ifdef SIFT3D_USE_OPENCL
	if (im_load_cl(cur, SIFT3D_FALSE))
		return SIFT3D_FAILURE;	
#endif

	f = (Sep_FIR_filter *) &gss->first_gauss.f;
	if (apply_Sep_FIR_filter(prev, cur, f, unit))
		return SIFT3D_FAILURE;

	// Build the rest of the pyramid
	SIFT3D_PYR_LOOP_LIMITED_START(o, s, o_start, o_end, s_start, s_end)
			cur = SIFT3D_PYR_IM_GET(gpyr, o, s);
			prev = SIFT3D_PYR_IM_GET(gpyr, o, s - 1);
			f = &gss->gauss_octave[s].f;
			if (apply_Sep_FIR_filter(prev, cur, f, unit))
				return SIFT3D_FAILURE;
#ifdef SIFT3D_USE_OPENCL
			if (im_read_back(cur, SIFT3D_FALSE))
				return SIFT3D_FAILURE;
#endif
		SIFT3D_PYR_LOOP_SCALE_END
		// Downsample
		if (o != o_end) {

                        const int downsample_level = 
                                SIFT3D_MAX(s_end - 2, gpyr->first_level);

			prev = SIFT3D_PYR_IM_GET(gpyr, o, downsample_level);
			cur = SIFT3D_PYR_IM_GET(gpyr, o + 1, s_start - 1);

                        assert(fabs(prev->s - cur->s) < FLT_EPSILON);

			if (im_downsample_2x(prev, cur))
				return SIFT3D_FAILURE;

		}
	SIFT3D_PYR_LOOP_OCTAVE_END

#ifdef SIFT3D_USE_OPENCL
	clFinish_all();
#endif

	return SIFT3D_SUCCESS;
}

static int build_dog(SIFT3D *sift3d) {

	Image *gpyr_cur, *gpyr_next, *dog_level;
	int o, s;

	Pyramid *const dog = &sift3d->dog;
	Pyramid *const gpyr = &sift3d->gpyr;

	SIFT3D_PYR_LOOP_START(dog, o, s)
		gpyr_cur = SIFT3D_PYR_IM_GET(gpyr, o, s);
		gpyr_next = SIFT3D_PYR_IM_GET(gpyr, o, s + 1);			
		dog_level = SIFT3D_PYR_IM_GET(dog, o, s);
		
		if (im_subtract(gpyr_cur, gpyr_next, 
						dog_level))
			return SIFT3D_FAILURE;
	SIFT3D_PYR_LOOP_END

	return SIFT3D_SUCCESS;
}

/* Detect local extrema */
static int detect_extrema(SIFT3D *sift3d, Keypoint_store *kp) {

	Image *cur, *prev, *next;
	Keypoint *key;
	float pcur, dogmax, peak_thresh;
	int o, s, x, y, z, x_start, x_end, y_start, y_end, z_start,
		z_end, num;

	const Pyramid *const dog = &sift3d->dog;
	const int o_start = dog->first_octave;
	const int o_end = SIFT3D_PYR_LAST_OCTAVE(dog);
	const int s_start = dog->first_level + 1;
	const int s_end = SIFT3D_PYR_LAST_LEVEL(dog) - 1;

	// Verify the inputs
	if (dog->num_levels < 3) {
		printf("detect_extrema: Requires at least 3 levels per octave, "
			   "provided only %d \n", dog->num_levels);
		return SIFT3D_FAILURE;
	}

	// Initialize dimensions of keypoint store
	cur = SIFT3D_PYR_IM_GET(dog, o_start, s_start);
	kp->nx = cur->nx;
	kp->ny = cur->ny;
	kp->nz = cur->nz;

#define CMP_CUBE(im, x, y, z, CMP, IGNORESELF, val) ( \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y),     (z) - 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y),     (z) - 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y),     (z) - 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) - 1, (z) - 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) + 1, (z) - 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) - 1, (z) - 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) - 1, (z) - 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) + 1, (z) - 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) + 1, (z) - 1, 0) && \
	((val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y),    (z), 0   ) || \
	    IGNORESELF) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y),     (z), 0    ) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y),     (z), 0    ) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) - 1, (z), 0    ) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) + 1, (z), 0    ) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) - 1, (z), 0    ) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) - 1, (z), 0    ) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) + 1, (z), 0    ) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) + 1, (z), 0    ) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y),     (z) + 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y),     (z) + 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y),     (z) + 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) - 1, (z) + 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) + 1, (z) + 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) - 1, (z) + 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) - 1, (z) + 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y) + 1, (z) + 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y) + 1, (z) + 1, 0) )
#ifdef CUBOID_EXTREMA
#define CMP_PREV(im, x, y, z, CMP, val) \
        CMP_CUBE(im, x, y, z, CMP, SIFT3D_FALSE, val)
#define CMP_CUR(im, x, y, z, CMP, val) \
        CMP_CUBE(im, x, y, z, CMP, SIFT3D_TRUE, val)
#define CMP_NEXT(im, x, y, z, CMP, val) \
        CMP_CUBE(im, x, y, z, CMP, SIFT3D_FALSE, val)
#else
#define CMP_PREV(im, x, y, z, CMP, val) ( \
        (val) CMP SIFT3D_IM_GET_VOX( (im), (x), (y), (z), 0) \
)
#define CMP_CUR(im, x, y, z, CMP, val) ( \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) + 1, (y),     (z), 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x) - 1, (y),     (z), 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) + 1, (z), 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y) - 1, (z), 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y),     (z) - 1, 0) && \
	(val) CMP SIFT3D_IM_GET_VOX( (im), (x),     (y),     (z) + 1, 0) \
)
#define CMP_NEXT(im, x, y, z, CMP, val) \
        CMP_PREV(im, x, y, z, CMP, val)
#endif

	num = 0;
	SIFT3D_PYR_LOOP_LIMITED_START(o, s, o_start, o_end, s_start, s_end)  

		// Select current and neighboring levels
		prev = SIFT3D_PYR_IM_GET(dog, o, s - 1);
		cur = SIFT3D_PYR_IM_GET(dog, o, s);
		next = SIFT3D_PYR_IM_GET(dog, o, s + 1);

		// Find maximum DoG value at this level
		dogmax = 0.0f;
		SIFT3D_IM_LOOP_START(cur, x, y, z)
			dogmax = SIFT3D_MAX(dogmax, 
                                fabsf(SIFT3D_IM_GET_VOX(cur, x, y, z, 0)));
		SIFT3D_IM_LOOP_END

		// Adjust threshold
		peak_thresh = sift3d->peak_thresh * dogmax;

		// Loop through all non-boundary pixels
		x_start = y_start = z_start = 1;
		x_end = cur->nx - 2;
		y_end = cur->ny - 2;
		z_end = cur->nz - 2;
		SIFT3D_IM_LOOP_LIMITED_START(cur, x, y, z, x_start, x_end, y_start,
							  y_end, z_start, z_end)
			// Sample the center value
			pcur = SIFT3D_IM_GET_VOX(cur, x, y, z, 0);

			// Apply the peak threshold
			if ((pcur > peak_thresh || pcur < -peak_thresh) && ((
				// Compare to the neighbors
				CMP_PREV(prev, x, y, z, >, pcur) &&
				CMP_CUR(cur, x, y, z, >, pcur) &&
				CMP_NEXT(next, x, y, z, >, pcur)
				) || (
				CMP_PREV(prev, x, y, z, <, pcur) &&
				CMP_CUR(cur, x, y, z, <, pcur) &&
				CMP_NEXT(next, x, y, z, <, pcur))))
				{

                                // Add a keypoint candidate
                                num++;
                                if (resize_Keypoint_store(kp, num))
                                        return SIFT3D_FAILURE;
                                key = kp->buf + num - 1;
                                if (init_Keypoint(key))
                                        return SIFT3D_FAILURE;
                                key->o = o;
                                key->s = s;
                                key->sd = cur->s;
				key->xd = (double) x;
				key->yd = (double) y;
				key->zd = (double) z;
                        }
		SIFT3D_IM_LOOP_END
	SIFT3D_PYR_LOOP_END
#undef CMP_NEIGHBORS

	return SIFT3D_SUCCESS;
}

/* Bin a Cartesian gradient into Spherical gradient bins */
SIFT3D_IGNORE_UNUSED
static int Cvec_to_sbins(const Cvec * const vd, Svec * const bins) {

	// Convert to spherical coordinates
	SIFT3D_CVEC_TO_SVEC(vd, bins);
	//FIXME: Is this needed? SIFT3D_CVEC_TO_SVEC cannot divide by zero
	if (bins->mag < FLT_EPSILON * 1E2)
		return SIFT3D_FAILURE;

	// Compute bins
	bins->az *= (float) NBINS_AZ / SIFT3D_AZ_MAX_F; 
	bins->po *= (float) NBINS_PO / SIFT3D_PO_MAX_F;

	assert(bins->az < NBINS_AZ);
	assert(bins->po < NBINS_PO);

	return SIFT3D_SUCCESS;
}

/* Refine a gradient histogram with optional operations,
 * such as solid angle weighting. */
static void refine_Hist(Hist *hist) {

#ifndef ICOS_HIST

#ifdef SIFT3D_ORI_SOLID_ANGLE_WEIGHT
	{	
	float po;
	int a, p;
	// TODO: could accelerate this with a lookup table		

	// Weight by the solid angle of the bins, ignoring constants
	HIST_LOOP_START(a, p)
		po = p * po_max_f / NBINS_PO;
		HIST_GET(hist, a, p) /= cosf(po) - cosf(po + 
			po_max_f / NBINS_PO);
	HIST_LOOP_END
	}
#endif

#endif

}

/* Assign rotation matrices to the keypoints. 
 * 
 * Note that this stage will modify kp, likely
 * rejecting some keypoints as orientationally
 * unstable. */
static int assign_orientations(SIFT3D *const sift3d, 
			       Keypoint_store *const kp) {

	Keypoint *kp_pos;
	size_t num;
	int i, err; 

	// Iterate over the keypoints 
        err = SIFT3D_SUCCESS;
#pragma omp parallel for
	for (i = 0; i < kp->slab.num; i++) {

		Keypoint *const key = kp->buf + i;
		const Image *const level = 
                        SIFT3D_PYR_IM_GET(&sift3d->gpyr, key->o, key->s);
                Mat_rm *const R = &key->R;
                const Cvec vcenter = {key->xd, key->yd, key->zd};
                const double sigma = ori_sig_fctr * key->sd;

		// Compute dominant orientations
                assert(R->u.data_float == key->r_data);
		switch (assign_orientation_thresh(level, &vcenter, sigma,
                                       sift3d->corner_thresh, R)) {
			case SIFT3D_SUCCESS:
				// Continue processing this keypoint
				break;
			case REJECT:
				// Mark this keypoint as invalid
                                key->xd = key->yd = key->zd = -1.0;
                                continue;
			default:
				// Any other return value is an error
                                err = SIFT3D_FAILURE;
                                continue;
		}
		
	}

        // Check for errors
        if (err) return err;

        // Rebuild the keypoint buffer in place
	kp_pos = kp->buf;
        for (i = 0; i < kp->slab.num; i++) {

		Keypoint *const key = kp->buf + i;

                // Check if the keypoint is valid
                if (key->xd < 0.0)
                        continue;

                // Copy this keypoint to the next available spot
                if (copy_Keypoint(key, kp_pos))
                        return SIFT3D_FAILURE;
               
                kp_pos++;
        }

	// Release unneeded keypoint memory
	num = kp_pos - kp->buf;
        return resize_Keypoint_store(kp, num);
}

/* Helper function to call assign_eig_ori, and reject keypoints with
 * confidence below the parameter "thresh." All other parameters are the same.
 * All return values are the same, except REJECT is returned if 
 * conf < thresh. */
static int assign_orientation_thresh(const Image *const im, 
        const Cvec *const vcenter, const double sigma, const double thresh,
        Mat_rm *const R) {

        double conf;
        int ret;

        ret = assign_eig_ori(im, vcenter, sigma, R, &conf);

        return ret == SIFT3D_SUCCESS ? 
                (conf < thresh ? REJECT : SIFT3D_SUCCESS) : ret;
}


/* Assign an orientation to a point in an image.
 *
 * Parameters:
 *   -im: The image data.
 *   -vcenter: The center of the window, in image space.
 *   -sigma: The scale parameter. The width of the window is a constant
 *      multiple of this.
 *   -R: The place to write the rotation matrix.
 */
static int assign_eig_ori(const Image *const im, const Cvec *const vcenter,
                          const double sigma, Mat_rm *const R, 
                          double *const conf) {

    Cvec v[2];
    Mat_rm A, L, Q;
    Cvec vd_win, vdisp, vr;
    double d, cos_ang, abs_cos_ang, corner_score;
    float weight, sq_dist, sgn;
    int i, x, y, z, m;
  
    const double win_radius = sigma * ori_rad_fctr; 

    // Verify inputs
    if (!SIFT3D_IM_CONTAINS_CVEC(im, vcenter)) {
        SIFT3D_ERR("assign_eig_ori: vcenter (%f, %f, %f) lies "
                "outside the boundaries of im [%d x %d x %d] \n", 
                vcenter->x, vcenter->y, vcenter->z, im->nx, im->ny, im->nz);
        return SIFT3D_FAILURE;
    }
    if (sigma < 0) {
        SIFT3D_ERR("assign_eig_ori: invalid sigma: %f \n", sigma);
        return SIFT3D_FAILURE;
    }

    // Initialize the intermediates
    if (init_Mat_rm(&A, 3, 3, SIFT3D_DOUBLE, SIFT3D_TRUE))
        return SIFT3D_FAILURE;
    if (init_Mat_rm(&L, 0, 0, SIFT3D_DOUBLE, SIFT3D_TRUE) ||
	init_Mat_rm(&Q, 0, 0, SIFT3D_DOUBLE, SIFT3D_TRUE))
	goto eig_ori_fail;

    // Resize the output
    R->num_rows = R->num_cols = IM_NDIMS;
    R->type = SIFT3D_FLOAT;
    if (resize_Mat_rm(R))
        goto eig_ori_fail;

    // Form the structure tensor and window gradient
    vd_win.x = 0.0f;
    vd_win.y = 0.0f;
    vd_win.z = 0.0f;
    IM_LOOP_SPHERE_START(im, x, y, z, vcenter, win_radius, &vdisp, sq_dist)

        Cvec vd;

	// Compute Gaussian weighting, ignoring the constant factor
	weight = expf(-0.5 * sq_dist / (sigma * sigma));		

	// Get the gradient	
	IM_GET_GRAD_ISO(im, x, y, z, 0, &vd);

	// Update the structure tensor
	SIFT3D_MAT_RM_GET(&A, 0, 0, double) += (double) vd.x * vd.x * weight;
	SIFT3D_MAT_RM_GET(&A, 0, 1, double) += (double) vd.x * vd.y * weight;
	SIFT3D_MAT_RM_GET(&A, 0, 2, double) += (double) vd.x * vd.z * weight;
	SIFT3D_MAT_RM_GET(&A, 1, 1, double) += (double) vd.y * vd.y * weight;
	SIFT3D_MAT_RM_GET(&A, 1, 2, double) += (double) vd.y * vd.z * weight;
	SIFT3D_MAT_RM_GET(&A, 2, 2, double) += (double) vd.z * vd.z * weight;

	// Update the window gradient
        SIFT3D_CVEC_SCALE(&vd, weight);
	SIFT3D_CVEC_OP(&vd_win, &vd, +, &vd_win);

    IM_LOOP_SPHERE_END

    // Fill in the remaining elements
    SIFT3D_MAT_RM_GET(&A, 1, 0, double) = SIFT3D_MAT_RM_GET(&A, 0, 1, double);
    SIFT3D_MAT_RM_GET(&A, 2, 0, double) = SIFT3D_MAT_RM_GET(&A, 0, 2, double);
    SIFT3D_MAT_RM_GET(&A, 2, 1, double) = SIFT3D_MAT_RM_GET(&A, 1, 2, double);

    // Reject keypoints with weak gradient 
    if (SIFT3D_CVEC_L2_NORM_SQ(&vd_win) < (float) ori_grad_thresh) {
	goto eig_ori_reject;
    } 

    // Get the eigendecomposition
    if (eigen_Mat_rm(&A, &Q, &L))
	goto eig_ori_fail;

    // Ensure we have distinct eigenvalues
    m = L.num_rows;
    if (m != 3)
	goto eig_ori_reject;

    // Test the eigenvectors for stability
    for (i = 0; i < m - 1; i++) {
	if (fabs(SIFT3D_MAT_RM_GET(&L, i, 0, double) /
		 SIFT3D_MAT_RM_GET(&L, i + 1, 0, double)) > max_eig_ratio)
	    goto eig_ori_reject;
    }

    // Assign signs to the first n - 1 vectors
    corner_score = DBL_MAX;
    for (i = 0; i < m - 1; i++) {

	const int eig_idx = m - i - 1;

	// Get an eigenvector, in descending order
	vr.x = (float) SIFT3D_MAT_RM_GET(&Q, 0, eig_idx, double);
	vr.y = (float) SIFT3D_MAT_RM_GET(&Q, 1, eig_idx, double);
	vr.z = (float) SIFT3D_MAT_RM_GET(&Q, 2, eig_idx, double);

	// Get the directional derivative
	d = SIFT3D_CVEC_DOT(&vd_win, &vr);

        // Get the cosine of the angle between the eigenvector and the gradient
        cos_ang = d / (SIFT3D_CVEC_L2_NORM(&vr) * SIFT3D_CVEC_L2_NORM(&vd_win));
        abs_cos_ang = fabs(cos_ang);

        // Compute the corner confidence score
        corner_score = SIFT3D_MIN(corner_score, abs_cos_ang);

	// Get the sign of the derivative
        sgn = d > 0.0 ? 1.0f : -1.0f;

	// Enforce positive directional derivative
	SIFT3D_CVEC_SCALE(&vr, sgn);

	// Add the vector to the rotation matrix
	SIFT3D_MAT_RM_GET(R, 0, i, float) = vr.x;
	SIFT3D_MAT_RM_GET(R, 1, i, float) = vr.y;
	SIFT3D_MAT_RM_GET(R, 2, i, float) = vr.z;

	// Save this vector for later use
	v[i] = vr;
    }

    // Take the cross product of the first two vectors
    SIFT3D_CVEC_CROSS(v, v + 1, &vr);

    // Add the last vector
    SIFT3D_MAT_RM_GET(R, 0, 2, float) = (float) vr.x;
    SIFT3D_MAT_RM_GET(R, 1, 2, float) = (float) vr.y;
    SIFT3D_MAT_RM_GET(R, 2, 2, float) = (float) vr.z;

    // Optionally write back the corner score
    if (conf != NULL)
        *conf = corner_score;

    cleanup_Mat_rm(&A);
    cleanup_Mat_rm(&Q);
    cleanup_Mat_rm(&L);
    return SIFT3D_SUCCESS; 

eig_ori_reject:
    if (conf != NULL)
        *conf = 0.0;
    cleanup_Mat_rm(&A);
    cleanup_Mat_rm(&Q);
    cleanup_Mat_rm(&L);
    return REJECT;

eig_ori_fail:
    if (conf != NULL)
        *conf = 0.0;
    cleanup_Mat_rm(&A);
    cleanup_Mat_rm(&Q);
    cleanup_Mat_rm(&L);
    return SIFT3D_FAILURE;
}

/* Assign an orientation to a point in an image.
 *
 * Parameters:
 *   -im: The image data.
 *   -kp: A container holding the keypoints. On successful return, the matrices
 *      "R" in these keypoints will be set to the assigned orientations.
 *   -conf: If not NULL, this is a pointer to array where the confidence scores
 *      are written. The scores are normally in the interval [0, 1], where 1 
 *      is the most confident, 0 is the least. A higher score means the assigned
 *      orientation is more likely to be robust. A negative score means the
 *      orientation could not be assigned. 
 *
 * If not NULL, this function will re-allocate conf. As such, *conf must either
 * be NULL or a pointer to a previously-allocated array.
 *
 * Return value: 
 * SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise.
 */
int SIFT3D_assign_orientations(const SIFT3D *const sift3d, 
        const Image *const im, Keypoint_store *const kp, double **const conf) {

        Image im_smooth;
        Keypoint key_base;
        int i;

        const int num = kp->slab.num;

        // Verify inputs 
        if (verify_keys(kp, im))
                return SIFT3D_FAILURE;

        // Initialize intermediates
        init_im(&im_smooth);
        init_Keypoint(&key_base);

        // Resize conf (num cannot be zero)
        if ((*conf = SIFT3D_safe_realloc(*conf, num * sizeof(double))) == NULL)
                goto assign_orientations_quit;

        // Smooth and scale the input
        if (smooth_scale_raw_input(sift3d, im, &im_smooth))
                goto assign_orientations_quit;

        // Assign each orientation
        for (i = 0; i < num; i++) {

                Cvec vcenter;
                int j;

                Keypoint *const key = kp->buf + i;
                double *const conf_ret = *conf + i;
                Mat_rm *const R = &key->R;

                // Convert the keypoint to the base octave
                if (keypoint2base(key, &key_base))
                        goto assign_orientations_quit;

                // Convert the keypoint to a vector
                vcenter.x = key_base.xd;
                vcenter.y = key_base.yd;
                vcenter.z = key_base.zd;

                // Assign the orientation
                switch (assign_eig_ori(&im_smooth, &vcenter, key_base.sd, R, 
                                       conf_ret))
                {
                        case SIFT3D_SUCCESS:
                                break;
                        case REJECT:
                                // Set R to identity
                                if (identity_Mat_rm(IM_NDIMS, R))
                                        goto assign_orientations_quit;
                                *conf_ret = -1.0;
                                break;
                        default:
                                // Stop processing if an error occurs
                                goto assign_orientations_quit;
                }
        }

        // Clean up
        im_free(&im_smooth);

        return SIFT3D_SUCCESS;

assign_orientations_quit:
        im_free(&im_smooth);
        return SIFT3D_FAILURE;
}

/* Detect keypoint locations and orientations. You must initialize
 * the SIFT3D struct, image, and keypoint store with the appropriate
 * functions prior to calling this function. */
int SIFT3D_detect_keypoints(SIFT3D *const sift3d, const Image *const im,
			    Keypoint_store *const kp) {

        // Verify inputs
        if (im->nc != 1) {
                SIFT3D_ERR("SIFT3D_detect_keypoints: invalid number "
                        "of image channels: %d -- only single-channel images "
                        "are supported \n", im->nc);
                return SIFT3D_FAILURE;
        }

        // Set the image       
        if (set_im_SIFT3D(sift3d, im))
                return SIFT3D_FAILURE;

	// Build the GSS pyramid
	if (build_gpyr(sift3d))
		return SIFT3D_FAILURE;

	// Build the DoG pyramid
	if (build_dog(sift3d))
		return SIFT3D_FAILURE;

	// Detect extrema
	if (detect_extrema(sift3d, kp))
		return SIFT3D_FAILURE;

	// Assign orientations
	if (assign_orientations(sift3d, kp))
		return SIFT3D_FAILURE;

	return SIFT3D_SUCCESS;
}

/* Get the bin and barycentric coordinates of a vector in the icosahedral 
 * histogram. */
SIFT3D_IGNORE_UNUSED
static int icos_hist_bin(const SIFT3D * const sift3d,
			   const Cvec * const x, Cvec * const bary,
			   int * const bin) { 

	float k;
	int i;

	const Mesh * const mesh = &sift3d->mesh;

	// Check for very small vectors
	if (SIFT3D_CVEC_L2_NORM_SQ(x) < bary_eps)
		return SIFT3D_FAILURE;

	// Iterate through the faces
	for (i = 0; i < ICOS_NFACES; i++) {

		const Tri * const tri = mesh->tri + i;

		// Convert to barycentric coordinates
		if (cart2bary(x, tri, bary, &k))
			continue;

		// Test for intersection
		if (bary->x < -bary_eps || bary->y < -bary_eps ||
		    bary->z < -bary_eps || k < 0)
			continue;

		// Save the bin
		*bin = i;

		// No other triangles will be intersected
		return SIFT3D_SUCCESS;
	}	

	// Unreachable code
	assert(SIFT3D_FALSE);
	return SIFT3D_FAILURE;
}

/* Helper routine to interpolate over the histograms of a
 * SIFT3D descriptor. */
void SIFT3D_desc_acc_interp(const SIFT3D * const sift3d, 
				const Cvec * const vbins, 
				const Cvec * const grad,
				SIFT3D_Descriptor * const desc) {

	Cvec dvbins;
	Hist *hist;
	float weight;
	int dx, dy, dz, x, y, z;

#ifdef ICOS_HIST
	Cvec bary;
	float mag;
	int bin;	
#else
	Svec sbins, dsbins;
	int da, dp, a, p;
#endif

	const int y_stride = NHIST_PER_DIM;
	const int z_stride = NHIST_PER_DIM * NHIST_PER_DIM; 

	// Compute difference from integer bin values
	dvbins.x = vbins->x - floorf(vbins->x);
	dvbins.y = vbins->y - floorf(vbins->y);
	dvbins.z = vbins->z - floorf(vbins->z);

	// Compute the histogram bin
#ifdef ICOS_HIST
	const Mesh *const mesh = &sift3d->mesh;

	// Get the index of the intersecting face 
	if (icos_hist_bin(sift3d, grad, &bary, &bin))
		return;
	
	// Get the magnitude of the vector
	mag = SIFT3D_CVEC_L2_NORM(grad);

#else
	if (Cvec_to_sbins(grad, &sbins))
		return;
	dsbins.az = sbins.az - floorf(sbins.az);
	dsbins.po = sbins.po - floorf(sbins.po);
#endif
	
	for (dx = 0; dx < 2; dx++) {
	for (dy = 0; dy < 2; dy++) {
        for (dz = 0; dz < 2; dz++) {

                x = (int) vbins->x + dx;
                y = (int) vbins->y + dy;
                z = (int) vbins->z + dz;

                // Check boundaries
                if (x < 0 || x >= NHIST_PER_DIM ||
                        y < 0 || y >= NHIST_PER_DIM ||
                        z < 0 || z >= NHIST_PER_DIM)
                        continue;

                // Get the histogram
                hist = desc->hists + x + y * y_stride + 
                           z * z_stride;	

                assert(x + y * y_stride + z * z_stride < DESC_NUM_TOTAL_HIST);

                // Get the spatial interpolation weight
                weight = ((dx == 0) ? (1.0f - dvbins.x) : dvbins.x) *
                        ((dy == 0) ? (1.0f - dvbins.y) : dvbins.y) *
                        ((dz == 0) ? (1.0f - dvbins.z) : dvbins.z);

                /* Add the value into the histogram */
#ifdef ICOS_HIST
                assert(HIST_NUMEL == ICOS_NVERT);
                assert(bin >= 0 && bin < ICOS_NFACES);

                // Interpolate over three vertices
                MESH_HIST_GET(mesh, hist, bin, 0) += mag * weight * bary.x;
                MESH_HIST_GET(mesh, hist, bin, 1) += mag * weight * bary.y;
                MESH_HIST_GET(mesh, hist, bin, 2) += mag * weight * bary.z; 
#else
                // Iterate over all angles
                for (dp = 0; dp < 2; dp ++) {
                for (da = 0; da < 2; da ++) {

                        a = ((int) sbins.az + da) % NBINS_AZ;
                        p = (int) sbins.po + dp;
                        if (p >= NBINS_PO) {
                                // See HIST_GET_PO
                                a = (a + NBINS_AZ / 2) % NBINS_AZ;
                                p = NBINS_PO - 1;
                        }
		
                        assert(a >= 0);
                        assert(a < NBINS_AZ);
                        assert(p >= 0);
                        assert(p < NBINS_PO);

                        HIST_GET(hist, a, p) += sbins.mag * weight *
                                ((da == 0) ? (1.0f - dsbins.az) : dsbins.az) *
                                ((dp == 0) ? (1.0f - dsbins.po) : dsbins.po);
                }}
#endif
	}}}

}

/* Normalize a descriptor */
static void normalize_desc(SIFT3D_Descriptor * const desc) {

	double norm; 
	int i, a, p;

	norm = 0.0;
	for (i = 0; i < DESC_NUM_TOTAL_HIST; i++) { 

                const Hist *const hist = desc->hists + i;

		HIST_LOOP_START(a, p) 
			const float el = HIST_GET(hist, a, p);
			norm += (double) el * el;
		HIST_LOOP_END 
	}

	norm = sqrt(norm) + DBL_EPSILON; 

	for (i = 0; i < DESC_NUM_TOTAL_HIST; i++) {

                Hist *const hist = desc->hists + i;
		const float norm_inv = 1.0f / norm; 

		HIST_LOOP_START(a, p) 
			HIST_GET(hist, a, p) *= norm_inv; 
		HIST_LOOP_END 
	}
}

/* Set a histogram to zero */
static void hist_zero(Hist *hist) {

        int a, p;

        HIST_LOOP_START(a, p)
                HIST_GET(hist, a, p) = 0.0f;
        HIST_LOOP_END
}

/* Helper routine to extract a single SIFT3D descriptor */
static int extract_descrip(SIFT3D *const sift3d, const Image *const im,
	   const Keypoint *const key, SIFT3D_Descriptor *const desc) {

        float buf[IM_NDIMS * IM_NDIMS];
        Mat_rm Rt;
	Cvec vcenter, vim, vkp, vbins, grad, grad_rot;
	Hist *hist;
	float weight, sq_dist;
	int i, x, y, z, a, p;

	// Compute basic parameters 
        const float sigma = key->sd * desc_sig_fctr;
	const float win_radius = desc_rad_fctr * sigma;
	const float desc_half_width = win_radius / sqrt(2);
	const float desc_width = 2.0f * desc_half_width;
        const float desc_hist_width = desc_width / NHIST_PER_DIM;
	const float desc_bin_fctr = 1.0f / desc_hist_width;
	const double coord_factor = ldexp(1.0, key->o);

        // Invert the rotation matrix
        if (init_Mat_rm_p(&Rt, buf, IM_NDIMS, IM_NDIMS, SIFT3D_FLOAT, 
		SIFT3D_FALSE) ||
                transpose_Mat_rm(&key->R, &Rt))
                return SIFT3D_FAILURE;

	// Zero the descriptor
	for (i = 0; i < DESC_NUM_TOTAL_HIST; i++) {
		hist = desc->hists + i;
                hist_zero(hist);
	}

	// Iterate over a sphere window in real-world coordinates 
	vcenter.x = key->xd;
	vcenter.y = key->yd;
	vcenter.z = key->zd;
	IM_LOOP_SPHERE_START(im, x, y, z, &vcenter, win_radius, &vim, sq_dist)

		// Rotate to keypoint space
		SIFT3D_MUL_MAT_RM_CVEC(&Rt, &vim, &vkp);		

		// Compute spatial bins
		vbins.x = (vkp.x + desc_half_width) * desc_bin_fctr;
		vbins.y = (vkp.y + desc_half_width) * desc_bin_fctr;
		vbins.z = (vkp.z + desc_half_width) * desc_bin_fctr;

		// Reject points outside the rectangular descriptor 
		if (vbins.x < 0 || vbins.y < 0 || vbins.z < 0 ||
			vbins.x >= (float) NHIST_PER_DIM ||
			vbins.y >= (float) NHIST_PER_DIM ||
			vbins.z >= (float) NHIST_PER_DIM)
			continue;

		// Take the gradient
		IM_GET_GRAD_ISO(im, x, y, z, 0, &grad);

		// Apply a Gaussian window
		weight = expf(-0.5f * sq_dist / (sigma * sigma));
		SIFT3D_CVEC_SCALE(&grad, weight);

                // Rotate the gradient to keypoint space
		SIFT3D_MUL_MAT_RM_CVEC(&Rt, &grad, &grad_rot);

		// Finally, accumulate bins by 5x linear interpolation
		SIFT3D_desc_acc_interp(sift3d, &vbins, &grad_rot, desc);
	IM_LOOP_SPHERE_END

	// Histogram refinement steps
	for (i = 0; i < DESC_NUM_TOTAL_HIST; i++) {
		refine_Hist(&desc->hists[i]);
	}

	// Normalize the descriptor
	normalize_desc(desc);

	// Truncate
	for (i = 0; i < DESC_NUM_TOTAL_HIST; i++) {
		hist = desc->hists + i;
		HIST_LOOP_START(a, p)
			HIST_GET(hist, a, p) = SIFT3D_MIN(HIST_GET(hist, a, p), 
						   (float) trunc_thresh);
		HIST_LOOP_END
	}

	// Normalize again
	normalize_desc(desc);

	// Save the descriptor location in the original image
	// coordinates
	desc->xd = key->xd * coord_factor;
	desc->yd = key->yd * coord_factor;
	desc->zd = key->zd * coord_factor;
	desc->sd = key->sd;

        return SIFT3D_SUCCESS;
}

/* Check if the Gaussian scale-space pyramid in a SIFT3D struct is valid. This
 * shall return SIFT3D_TRUE if the struct was initialized, and 
 * SIFT3D_detect_keypoints has been successfully called on it since 
 * initialization. 
 *
 * Note: sift3d must be initialized before calling this function. */
int SIFT3D_have_gpyr(const SIFT3D *const sift3d) {

        const Pyramid *const gpyr = &sift3d->gpyr;

        return gpyr->levels != NULL && gpyr->num_levels != 0 && 
                gpyr->num_octaves != 0;
}

/* Helper function to scale keypoint coordinates.
 *
 * Parameters:
 *  -src: The source keypoints.
 *  -factors: An array of length IM_NDIMS specifying the scaling factors.
 *  -dst: The scaled keypoints.
 *
 * Return: SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
static int scale_Keypoint(const Keypoint *const src, 
        const double *const factors, Keypoint *const dst) {

        // Copy the source keypoint
        if (copy_Keypoint(src, dst)) {
                SIFT3D_ERR("scale_Keypoint: failed to convert keypoints \n");
                return SIFT3D_FAILURE;
        }
      
        // Scale the coordinates 
        dst->xd *= factors[0];
        dst->yd *= factors[1];
        dst->zd *= factors[2];

        return SIFT3D_SUCCESS;
}

/* Helper function to smooth and scale a "raw" input image, as if it were 
 * processed via SIFT3D_detect_keypoints.
 *
 * Parameters:
 *  -sift3d: Stores the parameters sigma_n and sigma0.
 *  -src: The input "raw" image.
 *  -dst: The output, smoothed image.
 *
 * Return: SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise. */
static int smooth_scale_raw_input(const SIFT3D *const sift3d, 
        const Image *const src, Image *const dst) {

        Gauss_filter gauss;

        const double sigma_n = sift3d->gpyr.sigma_n;
        const double sigma0 = sift3d->gpyr.sigma0;
        const double unit = 1.0;

        // Initialize the smoothing filter        
        if (init_Gauss_incremental_filter(&gauss, sigma_n, sigma0, IM_NDIMS))
                return SIFT3D_FAILURE;

        // Smooth the input
        if (apply_Sep_FIR_filter(src, dst, &gauss.f, unit))
                goto smooth_scale_raw_input_quit;

        // Scale the input to [-1, 1]
        im_scale(dst);

        // Clean up
        cleanup_Gauss_filter(&gauss);
        
        return SIFT3D_SUCCESS;

smooth_scale_raw_input_quit:
        cleanup_Gauss_filter(&gauss);
        return SIFT3D_FAILURE;
}

/* Extract SIFT3D descriptors from a list of keypoints. Uses the Gaussian
 * scale-space pyramid from the previous call to SIFT3D_detect_keypoints on
 * this SIFT3D struct. To extract from an image, see 
 * SIFT3D_extract_raw_descriptors. 
 *
 * Note: To check if SIFT3D_detect_keypoints has been called on this struct,
 * use SIFT3D_have_gpyr.
 *
 * Parameters:
 *  sift3d - (initialized) struct defining the algorithm parameters. Must have
 *      been used in some previous call to SIFT3D_detect_keypoints.
 *  kp - keypoint list populated by a feature detector 
 *  desc - (initialized) struct to hold the descriptors
 *
 * Return value:
 *  Returns SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise.
 */
int SIFT3D_extract_descriptors(SIFT3D *const sift3d, 
        const Keypoint_store *const kp, 
        SIFT3D_Descriptor_store *const desc) {

	// Verify inputs
	if (verify_keys(kp, &sift3d->im))
		return SIFT3D_FAILURE;

        // Check if a Gaussian scale-space pyramid is available for processing
        if (!SIFT3D_have_gpyr(sift3d)) {
                SIFT3D_ERR("SIFT3D_extract_descriptors: no Gaussian pyramid is "
                        "available. Make sure SIFT3D_detect_keypoints was "
                        "called prior to calling this function. \n");
                return SIFT3D_FAILURE;
        }

        // Extract features
        if (_SIFT3D_extract_descriptors(sift3d, &sift3d->gpyr, kp, desc))
                return SIFT3D_FAILURE;

        return SIFT3D_SUCCESS;
}

/* Verify that keypoints kp are valid in image im. Returns SIFT3D_SUCCESS if
 * valid, SIFT3D_FAILURE otherwise. */
static int verify_keys(const Keypoint_store *const kp, const Image *const im) {

        int i;

	const int num = kp->slab.num;

	// Check the number of keypoints
	if (num < 1) {
		SIFT3D_ERR("verify_keys: invalid number of keypoints: "
				"%d \n", num);
		return SIFT3D_FAILURE;
	}

	// Check each keypoint
        for (i = 0; i < num; i++) {

                const Keypoint *key = kp->buf + i;

                const double octave_factor = ldexp(1.0, key->o);

                if (key->xd < 0 ||
                        key->yd < 0 ||
                        key->zd < 0 ||
                        key->xd * octave_factor >= (double) im->nx || 
                        key->yd * octave_factor >= (double) im->ny || 
                        key->zd * octave_factor >= (double) im->nz) {
                        SIFT3D_ERR("verify_keys: keypoint %d (%f, %f, %f) "
                                "octave %d exceeds image dimensions "
                                "(%d, %d, %d) \n", i, key->xd, key->yd, key->zd,
                                key->o, im->nx, im->ny, im->nz);
                        return SIFT3D_FAILURE; 
                }

                if (key->sd <= 0) {
                        SIFT3D_ERR("verify_keys: keypoint %d has invalid "
                                "scale %f \n", i, key->sd);
                        return SIFT3D_FAILURE;
                }
        }

        return SIFT3D_SUCCESS;
}

/* Convert the keypoint src to the equivalent one at octave 0, stored in dst. */
static int keypoint2base(const Keypoint *const src, Keypoint *const dst) {

        double base_factors[IM_NDIMS];
        int j;

        const double octave_factor = ldexp(1.0, src->o);

        // Convert the factors to the base octave
        for (j = 0; j < IM_NDIMS; j++) {
                base_factors[j] = octave_factor;
        }

        // Scale the keypoint
        if (scale_Keypoint(src, base_factors, dst))
                return SIFT3D_FAILURE;

        // Assign the keypoint the base octave and scale level
        dst->o = 0;
        dst->s = 0;

        return SIFT3D_SUCCESS;
}

/* Extract SIFT3D descriptors from a list of keypoints and 
 * an image. To extract from a Gaussian scale-space pyramid, see
 * SIFT3D_extract_descriptors. 
 *
 * Parameters:
 *  sift3d - (initialized) struct defining the algorithm parameters
 *  im - Pointer to an Image struct. A copy of the image is smoothed from 
 *      sift3d->sigma_n to sift3d->sigma0 prior to descriptor extraction.
 *  kp - keypoint list populated by a feature detector 
 *  desc - (initialized) struct to hold the descriptors
 *
 * Return value:
 *  Returns  SIFT3D_SUCCESS on success, SIFT3D_FAILURE otherwise.
 */
int SIFT3D_extract_raw_descriptors(SIFT3D *const sift3d, 
        const Image *const im, const Keypoint_store *const kp, 
        SIFT3D_Descriptor_store *const desc) {

        Keypoint_store kp_base;
        Pyramid pyr;
        Image *level;
        int i;

        const int first_octave = 0;
        const int first_level = 0;
        const int num_kp_levels = 1;
        const int num_levels = 1;
        const int num_octaves = 1;
        const double sigma0 = sift3d->gpyr.sigma0;
        const double sigma_n = sift3d->gpyr.sigma_n;

        // Verify inputs
        if (verify_keys(kp, im))
                return SIFT3D_FAILURE;

        // Initialize intermediates
        init_Pyramid(&pyr);
        init_Keypoint_store(&kp_base);

        // Initialize the image pyramid
        set_scales_Pyramid(sigma0, sigma_n, &pyr);
        if (resize_Pyramid(im, first_level, num_kp_levels, num_levels,
                first_octave, num_octaves, &pyr))
                goto extract_raw_descriptors_quit;

        // Smooth the input image, storing the result in the pyramid
        level = SIFT3D_PYR_IM_GET(&pyr, first_octave, first_level); 
        if (smooth_scale_raw_input(sift3d, im, level))
                goto extract_raw_descriptors_quit;

        // Allocate a temporary copy of the keypoints
        if (resize_Keypoint_store(&kp_base, kp->slab.num))
                goto extract_raw_descriptors_quit;

        // Convert keypoints to the base scale and octave
        for (i = 0; i < kp->slab.num; i++) {

                const Keypoint *const src = kp->buf + i;
                Keypoint *const dst = kp_base.buf + i; 

                if (keypoint2base(src, dst))
                        goto extract_raw_descriptors_quit;
        }

        // Extract the descriptors
        if (_SIFT3D_extract_descriptors(sift3d, &pyr, &kp_base, desc))
                goto extract_raw_descriptors_quit;

        // Clean up
        cleanup_Pyramid(&pyr);
        cleanup_Keypoint_store(&kp_base);

        return SIFT3D_SUCCESS;

extract_raw_descriptors_quit:
        cleanup_Pyramid(&pyr);
        cleanup_Keypoint_store(&kp_base);
        return SIFT3D_FAILURE;
}

/* Helper funciton to extract SIFT3D descriptors from a list of keypoints and 
 * an image. Called by SIFT3D_extract_descriptors and 
 * SIFT3D_extract_raw_descriptors.
 *
 * parameters:
 *  sift3d - (initialized) struct defining the algorithm parameters
 *  gpyr - A Gaussian Scale-Space pyramid containing the image data
 *  kp - keypoint list populated by a feature detector 
 *  desc - (initialized) struct to hold the descriptors
 *  use_gpyr - see im for details */
static int _SIFT3D_extract_descriptors(SIFT3D *const sift3d, 
        const Pyramid *const gpyr, const Keypoint_store *const kp, 
        SIFT3D_Descriptor_store *const desc) {

	int i, ret;

	const Image *const first_level = 
                SIFT3D_PYR_IM_GET(gpyr, gpyr->first_octave, gpyr->first_level);

	const int num = kp->slab.num;

	// Initialize the metadata 
	desc->nx = first_level->nx;	
	desc->ny = first_level->ny;	
	desc->nz = first_level->nz;	

	// Resize the descriptor store
        if (resize_SIFT3D_Descriptor_store(desc, num))
                return SIFT3D_FAILURE;

        // Extract the descriptors
        ret = SIFT3D_SUCCESS;
#pragma omp parallel for
	for (i = 0; i < desc->num; i++) {

                const Keypoint *const key = kp->buf + i;
		SIFT3D_Descriptor *const descrip = desc->buf + i;
		const Image *const level = 
                        SIFT3D_PYR_IM_GET(gpyr, key->o, key->s);

		if (extract_descrip(sift3d, level, key, descrip)) {
                        ret = SIFT3D_FAILURE;
                }
	}	

	return ret;
}

/* L2-normalize a histogram */
static void normalize_hist(Hist *const hist) {

        double norm;
        float norm_inv;
        int a, p;

        norm = 0.0;
        HIST_LOOP_START(a, p)
                const float el = HIST_GET(hist, a, p);
                norm += (double) el * el;
        HIST_LOOP_END

        norm = sqrt(norm) + DBL_EPSILON;
        norm_inv = 1.0f / norm; 

        HIST_LOOP_START(a, p)
                 HIST_GET(hist, a, p) *= norm_inv; 
        HIST_LOOP_END 
}

/* Dense gradient histogram postprocessing steps */
static void postproc_Hist(Hist *const hist, const float norm) {

        int a, p;

        const float hist_trunc = trunc_thresh * DESC_NUMEL / HIST_NUMEL;

	// Histogram refinement steps
	refine_Hist(hist);

	// Normalize the descriptor
	normalize_hist(hist);

	// Truncate
	HIST_LOOP_START(a, p)
		HIST_GET(hist, a, p) = SIFT3D_MIN(HIST_GET(hist, a, p), 
					   hist_trunc);
	HIST_LOOP_END

	// Normalize again
	normalize_hist(hist);

        // Convert to the desired norm
        HIST_LOOP_START(a, p)
                HIST_GET(hist, a, p) *= norm;
        HIST_LOOP_END
}

/* Helper routine to extract a single SIFT3D histogram, with rotation. */
static int extract_dense_descrip_rotate(SIFT3D *const sift3d, 
           const Image *const im, const Cvec *const vcenter, 
           const double sigma, const Mat_rm *const R, Hist *const hist) {

        float buf[IM_NDIMS * IM_NDIMS];
        Mat_rm Rt;
	Cvec grad, grad_rot, bary, vim;
	float sq_dist, mag, weight;
SIFT3D_IGNORE_UNUSED
	int a, p, x, y, z, bin;

        const Mesh *const mesh = &sift3d->mesh;
	const float win_radius = desc_rad_fctr * sigma;

        // Invert the rotation matrix
        if (init_Mat_rm_p(&Rt, buf, IM_NDIMS, IM_NDIMS, SIFT3D_FLOAT, 
		SIFT3D_FALSE) ||
                transpose_Mat_rm(R, &Rt))
                return SIFT3D_FAILURE;

	// Zero the descriptor
        hist_zero(hist);

	// Iterate over a sphere window in real-world coordinates
	IM_LOOP_SPHERE_START(im, x, y, z, vcenter, win_radius, &vim, sq_dist)

		// Take the gradient and rotate
		IM_GET_GRAD_ISO(im, x, y, z, 0, &grad);
		SIFT3D_MUL_MAT_RM_CVEC(&Rt, &grad, &grad_rot);

                // Get the index of the intersecting face
                if (icos_hist_bin(sift3d, &grad_rot, &bary, &bin))
                        continue;

                // Get the magnitude of the vector
                mag = SIFT3D_CVEC_L2_NORM(&grad);

		// Get the Gaussian window weight
		weight = expf(-0.5f * sq_dist / (sigma * sigma));

                // Interpolate over three vertices
                MESH_HIST_GET(mesh, hist, bin, 0) += mag * weight * bary.x;
                MESH_HIST_GET(mesh, hist, bin, 1) += mag * weight * bary.y;
                MESH_HIST_GET(mesh, hist, bin, 2) += mag * weight * bary.z;

	IM_LOOP_SPHERE_END

        return SIFT3D_SUCCESS;
}

/* Get a descriptor with a single histogram at each voxel of an image.
 * The result is an image with HIST_NUMEL channels, where each channel is a
 * bin of the histogram.
 *
 * Parameters:
 * -sift3d Stores the algorithm parameters.
 * -in The input image.
 * -desc The output image of descriptors.
 */
int SIFT3D_extract_dense_descriptors(SIFT3D *const sift3d, 
        const Image *const in, Image *const desc) {

        int (*extract_fun)(SIFT3D *const, const Image *const, Image *const);
        Image in_smooth;
        int x, y, z;

        // Verify inputs
        if (in->nc != 1) {
                SIFT3D_ERR("SIFT3D_extract_dense_descriptors: invalid "
                        "number of channels: %d. This function only supports "
                        "single-channel images. \n", in->nc);
                return SIFT3D_FAILURE;
        }

        // Select the appropriate subroutine
        extract_fun = sift3d->dense_rotate ? 
                extract_dense_descriptors_rotate :
                extract_dense_descriptors_no_rotate;

        // Resize the output image
        memcpy(SIFT3D_IM_GET_DIMS(desc), SIFT3D_IM_GET_DIMS(in), 
                IM_NDIMS * sizeof(int));
        desc->nc = HIST_NUMEL;
        im_default_stride(desc);
        if (im_resize(desc))
                return SIFT3D_FAILURE;

        // Intialize intermediates
        init_im(&in_smooth);

        //TODO: Interpolate to be isotropic

        // Smooth and scale the input image
        if (smooth_scale_raw_input(sift3d, in, &in_smooth))
                goto extract_dense_quit;

        // Extract the descriptors
        if (extract_fun(sift3d, &in_smooth, desc))
                return SIFT3D_FAILURE;

        // Post-process the descriptors
        SIFT3D_IM_LOOP_START(desc, x, y, z)

                Hist hist;

                // Get the image intensity at this voxel 
                const float val = SIFT3D_IM_GET_VOX(in, x, y, z, 0);

                // Copy to a Hist
                vox2hist(desc, x, y, z, &hist);

                // Post-process
                postproc_Hist(&hist, val);

                // Copy back to the image
                hist2vox(&hist, desc, x, y, z);

        SIFT3D_IM_LOOP_END

        // TODO transform back to original space

        // Clean up
        im_free(&in_smooth);

        return SIFT3D_SUCCESS;

extract_dense_quit:
        im_free(&in_smooth);
        return SIFT3D_FAILURE;
}

/* Helper function for SIFT3D_extract_dense_descriptors, without rotation 
 * invariance. This function is much faster than its rotation-invariant 
 * counterpart because histogram bins are pre-computed. */
static int extract_dense_descriptors_no_rotate(SIFT3D *const sift3d,
        const Image *const in, Image *const desc) {

        Image temp; 
        Gauss_filter gauss;
	Cvec grad, bary;
        int x, y, z, bin;

        const int x_start = 1;
        const int y_start = 1;
        const int z_start = 1;
        const int x_end = in->nx - 2;
        const int y_end = in->ny - 2;
        const int z_end = in->nz - 2;

        Mesh * const mesh = &sift3d->mesh;
        const double sigma_win = sift3d->gpyr.sigma0 * desc_sig_fctr / 
                                 NHIST_PER_DIM;
        const double unit = 1.0;

        // Initialize the intermediate image
        init_im(&temp);
        if (im_copy_dims(desc, &temp))
                return SIFT3D_FAILURE;

        // Initialize the filter
        if (init_Gauss_filter(&gauss, sigma_win, 3)) {
                im_free(&temp);
                return SIFT3D_FAILURE;
        }

        // Initialize the descriptors for each voxel
        im_zero(&temp);
        SIFT3D_IM_LOOP_LIMITED_START(in, x, y, z, x_start, x_end, y_start, 
                y_end, z_start, z_end)

                // Take the gradient
		IM_GET_GRAD_ISO(in, x, y, z, 0, &grad);

                // Get the index of the intersecting face
                if (icos_hist_bin(sift3d, &grad, &bary, &bin))
                        continue;

                // Initialize each vertex
                SIFT3D_IM_GET_VOX(&temp, x, y, z, MESH_GET_IDX(mesh, bin, 0)) = 
                        bary.x;
                SIFT3D_IM_GET_VOX(&temp, x, y, z, MESH_GET_IDX(mesh, bin, 1)) = 
                        bary.y;
                SIFT3D_IM_GET_VOX(&temp, x, y, z, MESH_GET_IDX(mesh, bin, 2)) = 
                        bary.z;

        SIFT3D_IM_LOOP_END

        // Filter the descriptors
	if (apply_Sep_FIR_filter(&temp, desc, &gauss.f, unit))
                goto dense_extract_quit;

        // Clean up
        im_free(&temp);
        cleanup_Gauss_filter(&gauss);

        return SIFT3D_SUCCESS;

dense_extract_quit:
        im_free(&temp);
        cleanup_Gauss_filter(&gauss);
        return SIFT3D_FAILURE;
}

/* Copy a voxel to a Hist. Does no bounds checking. */
static void vox2hist(const Image *const im, const int x, const int y,
        const int z, Hist *const hist) {

        int c;

        for (c = 0; c < HIST_NUMEL; c++) {
                hist->bins[c] = SIFT3D_IM_GET_VOX(im, x, y, z, c);
        }
}

/* Copy a Hist to a voxel. Does no bounds checking. */
static void hist2vox(Hist *const hist, const Image *const im, const int x, 
        const int y, const int z) {

        int c;
        
        for (c = 0; c < HIST_NUMEL; c++) {
                SIFT3D_IM_GET_VOX(im, x, y, z, c) = hist->bins[c];
        }
}

/* As in extract_dense_descrip, but with rotation invariance */
static int extract_dense_descriptors_rotate(SIFT3D *const sift3d,
        const Image *const in, Image *const desc) {

        Hist hist;
        Mat_rm R, Id;
        Mat_rm *ori;
        int i, x, y, z;

        // Initialize the identity matrix
        if (init_Mat_rm(&Id, 3, 3, SIFT3D_FLOAT, SIFT3D_TRUE)) {
                return SIFT3D_FAILURE;
        }
        for (i = 0; i < 3; i++) {       
                SIFT3D_MAT_RM_GET(&Id, i, i, float) = 1.0f;
        }

        // Initialize the rotation matrix
        if (init_Mat_rm(&R, 3, 3, SIFT3D_FLOAT, SIFT3D_TRUE)) {
                cleanup_Mat_rm(&Id);
                return SIFT3D_FAILURE;
        }

        // Iterate over each voxel
        SIFT3D_IM_LOOP_START(in, x, y, z)

                const Cvec vcenter = {x, y, z}; 

                const double ori_sigma = sift3d->gpyr.sigma0 * ori_sig_fctr;
                const double desc_sigma = sift3d->gpyr.sigma0 * 
                        desc_sig_fctr / NHIST_PER_DIM;

                // Attempt to assign an orientation
                switch (assign_orientation_thresh(in, &vcenter, ori_sigma,
                                      sift3d->corner_thresh, &R)) {
                        case SIFT3D_SUCCESS:
                                // Use the assigned orientation
                                ori = &R;
                                break;
                        case REJECT:
                                // Default to identity
                                ori = &Id;
                                break;
                        default:
                                // Unexpected error
                                goto dense_rotate_quit;
                }

                // Extract the descriptor
                if (extract_dense_descrip_rotate(sift3d, in, &vcenter, 
                        desc_sigma, ori, &hist))
                        goto dense_rotate_quit;

                // Copy the descriptor to the image channels
                hist2vox(&hist, desc, x, y, z);

        SIFT3D_IM_LOOP_END

        // Clean up
        cleanup_Mat_rm(&R);
        cleanup_Mat_rm(&Id);
        return SIFT3D_SUCCESS;

dense_rotate_quit:
        // Clean up and return an error condition 
        cleanup_Mat_rm(&R);
        cleanup_Mat_rm(&Id);
        return SIFT3D_FAILURE;
}

/* Convert a keypoint store to a matrix. 
 * Output format:
 *  [x1 y1 z1]
 *  |   ...  |
 *  [xn yn zn] 
 * 
 * mat must be initialized. */
int Keypoint_store_to_Mat_rm(const Keypoint_store *const kp, Mat_rm *const mat) {

    int i;

    const int num = kp->slab.num;

    // Resize mat
    mat->num_rows = num;
    mat->num_cols = IM_NDIMS;
    mat->type = SIFT3D_DOUBLE;
    if (resize_Mat_rm(mat))
	return SIFT3D_FAILURE;

    // Build the matrix
    for (i = 0; i < num; i++) {

        const Keypoint *const key = kp->buf + i;

        // Adjust the coordinates to the base octave
        const double coord_factor = ldexp(1.0, key->o);

	SIFT3D_MAT_RM_GET(mat, i, 0, double) = coord_factor * key->xd;
	SIFT3D_MAT_RM_GET(mat, i, 1, double) = coord_factor * key->yd;
	SIFT3D_MAT_RM_GET(mat, i, 2, double) = coord_factor * key->zd;
    }

    return SIFT3D_SUCCESS;
}

/* Convert SIFT3D desriptors to a coordinate matrix. This function has the same
 * output format as Keypoint_store_to_Mat_rm */
int SIFT3D_Descriptor_coords_to_Mat_rm(
        const SIFT3D_Descriptor_store *const store, 
        Mat_rm *const mat) {

        int i;

	const int num_rows = store->num;
        const int num_cols = IM_NDIMS;

	// Verify inputs
	if (num_rows < 1) {
		printf("SIFT3D_Descriptor_coords_to_Mat_rm: invalid number of "
		       "descriptors: %d \n", num_rows);
		return SIFT3D_FAILURE;
	}

        // Resize the output
        mat->type = SIFT3D_DOUBLE;
        mat->num_rows = num_rows;
        mat->num_cols = num_cols;
        if (resize_Mat_rm(mat))
                return SIFT3D_FAILURE;

        // Copy the data
        for (i = 0; i < num_rows; i++) {

		const SIFT3D_Descriptor *const desc = store->buf + i;

                SIFT3D_MAT_RM_GET(mat, i, 0, double) = desc->xd;
                SIFT3D_MAT_RM_GET(mat, i, 1, double) = desc->yd;
                SIFT3D_MAT_RM_GET(mat, i, 2, double) = desc->zd;
        }

        return SIFT3D_SUCCESS;
}

/* Convert SIFT3D descriptors to a matrix.
 *
 * Output format:
 *  [x y z el0 el1 ... el767]
 * Each row is a feature descriptor. [x y z] are the coordinates in image
 * space, and [el0 el1 ... el767 are the 768 dimensions of the descriptor.
 *
 * mat must be initialized prior to calling this function. mat will be resized.
 * The type of mat will be changed to float.
 */
int SIFT3D_Descriptor_store_to_Mat_rm(const SIFT3D_Descriptor_store *const store, 
				      Mat_rm *const mat) {
	int i, j, a, p;

	const int num_rows = store->num;
	const int num_cols = IM_NDIMS + DESC_NUMEL;

	// Verify inputs
	if (num_rows < 1) {
		printf("SIFT3D_Descriptor_store_to_Mat_rm: invalid number of "
		       "descriptors: %d \n", num_rows);
		return SIFT3D_FAILURE;
	}

	// Resize inputs
	mat->type = SIFT3D_FLOAT;
	mat->num_rows = num_rows;
	mat->num_cols = num_cols;
	if (resize_Mat_rm(mat))
		return SIFT3D_FAILURE;

	// Copy the data
	for (i = 0; i < num_rows; i++) {

		const SIFT3D_Descriptor *const desc = store->buf + i;

		// Copy the coordinates
		SIFT3D_MAT_RM_GET(mat, i, 0, float) = (float) desc->xd;
		SIFT3D_MAT_RM_GET(mat, i, 1, float) = (float) desc->yd;
		SIFT3D_MAT_RM_GET(mat, i, 2, float) = (float) desc->zd;

		// Copy the feature vector
		for (j = 0; j < DESC_NUM_TOTAL_HIST; j++) {
			const Hist *const hist = desc->hists + j;
			HIST_LOOP_START(a, p)
                                const int col = DESC_MAT_GET_COL(j, a, p);
				SIFT3D_MAT_RM_GET(mat, i, col, float) = 
					HIST_GET(hist, a, p);
			HIST_LOOP_END
		}
	}

	return SIFT3D_SUCCESS;
}

/* Convert a Mat_rm to a descriptor store. See 
 * SIFT3D_Descriptor_store_to_Mat_rm for the input format. */
int Mat_rm_to_SIFT3D_Descriptor_store(const Mat_rm *const mat, 
				      SIFT3D_Descriptor_store *const store) {

	int i, j, a, p;

	const int num_rows = mat->num_rows;
	const int num_cols = mat->num_cols;

	// Verify inputs
	if (num_rows < 1 || num_cols != IM_NDIMS + DESC_NUMEL) {
		SIFT3D_ERR("Mat_rm_to_SIFT3D_Descriptor_store: invalid matrix "
		       "dimensions: [%d X %d] \n", num_rows, num_cols);
		return SIFT3D_FAILURE;
	}
        if (mat->type != SIFT3D_FLOAT) {
		SIFT3D_ERR("Mat_rm_to_SIFT3D_Descriptor_store: matrix must "
                        "have type SIFT3D_FLOAT");
		return SIFT3D_FAILURE;
        }

        /* Resize the descriptor store */
        if (resize_SIFT3D_Descriptor_store(store, num_rows))
                return SIFT3D_FAILURE;

	// Copy the data
	for (i = 0; i < num_rows; i++) {

		SIFT3D_Descriptor *const desc = store->buf + i;

		// Copy the coordinates
		desc->xd = SIFT3D_MAT_RM_GET(mat, i, 0, float);
		desc->yd = SIFT3D_MAT_RM_GET(mat, i, 1, float);
		desc->zd = SIFT3D_MAT_RM_GET(mat, i, 2, float);
                desc->sd = sigma0_default;

		// Copy the feature vector
		for (j = 0; j < DESC_NUM_TOTAL_HIST; j++) {
			Hist *const hist = desc->hists + j;
			HIST_LOOP_START(a, p)
                                const int col = DESC_MAT_GET_COL(j, a, p);
				HIST_GET(hist, a, p) = 
                                        SIFT3D_MAT_RM_GET(mat, i, col, float);
			HIST_LOOP_END
		}
	}
	
	return SIFT3D_SUCCESS;
}

/* Convert a list of matches to matrices of point coordinates.
 * Only valid matches will be included in the output matrices.
 *
 * The format of "matches" is specified in SIFT3D_nn_match.
 *
 * All matrices must be initialized prior to calling this function.
 *
 * Output format:
 *  m x 3 matrices [x11 y11 z11] [x21 y21 z21]
 * 		   |x12 y12 z12| |x22 y22 z22|
 *		        ...	      ...
 * 		   [x1N y1N z1N] [x2N y2N z2N] 
 *
 * Where points on corresponding rows are matches. */
int SIFT3D_matches_to_Mat_rm(SIFT3D_Descriptor_store *d1,
			     SIFT3D_Descriptor_store *d2,
			     const int *const matches,
			     Mat_rm *const match1, 
			     Mat_rm *const match2) {
    int i, num_matches;

    const int num = d1->num;

    // Resize matrices 
    match1->num_rows = match2->num_rows = d1->num;
    match1->num_cols = match2->num_cols = 3;
    match1->type = match2->type = SIFT3D_DOUBLE;
    if (resize_Mat_rm(match1) || resize_Mat_rm(match2))
	    return SIFT3D_FAILURE;

    // Populate the matrices
    num_matches = 0;
    for (i = 0; i < num; i++) {

	    const SIFT3D_Descriptor *const desc1 = d1->buf + i;
	    const SIFT3D_Descriptor *const desc2 = d2->buf + matches[i];

	    if (matches[i] == -1)
		    continue;

	    // Save the match
	    SIFT3D_MAT_RM_GET(match1, num_matches, 0, double) = desc1->xd; 
	    SIFT3D_MAT_RM_GET(match1, num_matches, 1, double) = desc1->yd; 
	    SIFT3D_MAT_RM_GET(match1, num_matches, 2, double) = desc1->zd; 
	    SIFT3D_MAT_RM_GET(match2, num_matches, 0, double) = desc2->xd; 
	    SIFT3D_MAT_RM_GET(match2, num_matches, 1, double) = desc2->yd; 
	    SIFT3D_MAT_RM_GET(match2, num_matches, 2, double) = desc2->zd; 
	    num_matches++;
    }

    // Release extra memory
    match1->num_rows = match2->num_rows = num_matches;
    if (resize_Mat_rm(match1) || resize_Mat_rm(match2))
	    return SIFT3D_FAILURE;
    
    return SIFT3D_SUCCESS;
}

/* Perform nearest neighbor matching on two sets of 
 * SIFT descriptors.
 *
 * This function will reallocate *matches. As such, *matches must be either
 * NULL or a pointer to previously-allocated array. Upon successful exit,
 * *matches is an array of size d1->num.
 * 
 * On return, the ith element of matches contains the index in d2 of the match
 * corresponding to the ith descriptor in d1, or -1 if no match was found.
 *
 * You might consider using SIFT3D_matches_to_Mat_rm to convert the matches to
 * coordinate matrices. */
int SIFT3D_nn_match(const SIFT3D_Descriptor_store *const d1,
		    const SIFT3D_Descriptor_store *const d2,
		    const float nn_thresh, int **const matches) {

	int i;

	const int num = d1->num;

        // Verify inputs
	if (num < 1) {
		SIFT3D_ERR("_SIFT3D_nn_match: invalid number of "
			"descriptors in d1: %d \n", num);
		return SIFT3D_FAILURE;
	}

	// Resize the matches array (num cannot be zero)
	if ((*matches = (int *) SIFT3D_safe_realloc(*matches, 
		num * sizeof(int))) == NULL) {
	    SIFT3D_ERR("_SIFT3D_nn_match: out of memory! \n");
	    return SIFT3D_FAILURE;
	}

	for (i = 0; i < d1->num; i++) {
	    // Mark -1 to signal there is no match
	    (*matches)[i] = -1;
	}
	
	// Exhaustive search for matches
#pragma omp parallel for
	for (i = 0; i < num; i++) {

                const SIFT3D_Descriptor *const desc1 = d1->buf + i;
                int *const match = *matches + i;

                // Forward matching pass
                *match = match_desc(desc1, d2, nn_thresh);

                // We are done if there was no match
                if (*match < 0)
                        continue;

                // Check for forward-backward consistency
                if (match_desc(d2->buf + *match, d1, nn_thresh) != i) {
                        *match = -1;
                }
        }

	return SIFT3D_SUCCESS;
}

/* Helper function to match desc against the descriptors in store. Returns the
 * index of the match, or -1 if none was found. */
static int match_desc(const SIFT3D_Descriptor *const desc,
        const SIFT3D_Descriptor_store *const store, const float nn_thresh) {

	const SIFT3D_Descriptor *desc_best;
        double ssd_best, ssd_nearest;
        int i;

#ifdef SIFT3D_MATCH_MAX_DIST
        Cvec dims, dmatch;
        double dist_match;
				
        // Compute spatial distance rejection threshold
        dims.x = (float) store->nx;	
        dims.y = (float) store->ny;	
        dims.z = (float) store->nz;	
        const double diag = SIFT3D_CVEC_L2_NORM(&dims);	
        const double dist_thresh = diag * SIFT3D_MATCH_MAX_DIST;
#endif

        // Linear search for the best and second-best SSD matches 
        ssd_best = ssd_nearest = DBL_MAX;
        desc_best = NULL;
        for (i = 0; i < store->num; i++) { 

                double ssd;
                int j;

                const SIFT3D_Descriptor *const desc2 = store->buf + i;

                // Compute the SSD of the two descriptors
                ssd = 0.0;
                for (j = 0; j < DESC_NUM_TOTAL_HIST; j++) {

                        int a, p;

                        const Hist *const hist1 = desc->hists + j;
                        const Hist *const hist2 = desc2->hists + j;

                        HIST_LOOP_START(a, p)
                                const double diff = 
                                        (double) HIST_GET(hist1, a, p) -
                                        (double) HIST_GET(hist2, a, p);
                                        ssd += diff * diff;
                        HIST_LOOP_END

                        // Early termination
                        if (ssd > ssd_nearest)
                                break;
                }

                // Compare to the best matches
                if (ssd < ssd_best) {
                        desc_best = desc2; 
                        ssd_nearest = ssd_best;
                        ssd_best = ssd;
                } else  {
                        ssd_nearest = SIFT3D_MIN(ssd_nearest, ssd);
                }
        }

        // Reject a match if the nearest neighbor is too close
        if (ssd_best / ssd_nearest > nn_thresh * nn_thresh)
                        return -1;

#ifdef SIFT3D_MATCH_MAX_DIST
        // Compute the spatial distance of the match
        dmatch.x = (float) desc_best->xd - desc1->xd; 
        dmatch.y = (float) desc_best->yd - desc1->yd; 
        dmatch.z = (float) desc_best->zd - desc1->zd; 
        dist_match = (double) SIFT3D_CVEC_L2_NORM(&dmatch);

        // Reject matches of great distance
        if (dist_match > dist_thresh)
                return -1;
#endif
        // The match was a success
        return desc_best - store->buf;
}
			
/* Draw the matches. 
 * 
 * Inputs:
 * -left - lefthand-side image
 * -right - righthand-side image
 * -keys_left - Keypoints from "left" (can be NULL if keys is NULL) 
 * -keys_right - Keypoints from "right" (can be NULL if keys is NULL)
 * -match_left - Matches from "left" (can be NULL if lines is NULL)
 * -match_right - Matches from "right" (can be NULL if lines is NULL)
 * 
 * Outputs:
 * -concat - Concatenated image (NULL if not needed)
 * -keys - Keypoints from concat (NULL is not needed)
 * -lines - Lines between matching keypoints in concat (NULL if not needed)
 * It is an error if all outputs are NULL.
 *
 * Return:
 * SIFT3D_SUCCESS or SIFT3D_FAILURE
 */
int draw_matches(const Image *const left, const Image *const right,
                 const Mat_rm *const keys_left, const Mat_rm *const keys_right,
		 const Mat_rm *const match_left, 
                 const Mat_rm *const match_right, Image *const concat, 
                 Image *const keys, Image *const lines) {

        Image concat_temp, left_padded, right_padded;
	Mat_rm keys_right_draw, keys_left_draw, keys_draw, match_right_draw;
	int i;

        const double right_pad = (double) left->nx;
        const int ny_pad = SIFT3D_MAX(right->ny, left->ny);
	const int nz_pad = SIFT3D_MAX(right->nz, left->nz);

        // Choose which image to use for concatenation 
        Image *const concat_arg = concat == NULL ? &concat_temp : concat;

        // Verify inputs
        if (concat == NULL && keys == NULL && lines == NULL) {
                SIFT3D_ERR("draw_matches: all outputs are NULL \n");
                return SIFT3D_FAILURE;
        }
        if (keys_left == NULL && keys != NULL) {
                SIFT3D_ERR("draw_matches: keys_left is NULL but keys is "
                        "not \n");
                return SIFT3D_FAILURE;
        }
        if (keys_right == NULL && keys != NULL) {
                SIFT3D_ERR("draw_matches: keys_right is NULL but keys is "
                        "not \n");
                return SIFT3D_FAILURE;
        }
        if (match_left == NULL && lines != NULL) {
                SIFT3D_ERR("draw_matches: match_left is NULL but lines "
                        "is not \n");
                return SIFT3D_FAILURE;
        }
        if (match_right == NULL && lines != NULL) {
                SIFT3D_ERR("draw_matches: match_right is NULL but lines "
                        "is not \n");
                return SIFT3D_FAILURE;
        }

	// Initialize intermediates		
        init_im(&concat_temp);
        init_im(&left_padded);
        init_im(&right_padded);
        if (init_Mat_rm(&keys_right_draw, 0, 0, SIFT3D_DOUBLE, 
		SIFT3D_FALSE) ||
	        init_Mat_rm(&match_right_draw, 0, 0, SIFT3D_DOUBLE, 
			SIFT3D_FALSE) ||
                init_Mat_rm(&keys_left_draw, 0, 0, SIFT3D_DOUBLE, 
			SIFT3D_FALSE) ||
                init_Mat_rm(&keys_draw, 0, 0, SIFT3D_DOUBLE, SIFT3D_FALSE))
	        return SIFT3D_FAILURE;

        // Pad the images to be the same in all dimensions but x
	if (init_im_with_dims(&right_padded, right->nx, ny_pad, nz_pad, 1) || 
	        init_im_with_dims(&left_padded, left->nx, ny_pad, nz_pad, 1) ||
	   	im_pad(right, &right_padded) || 
	    	im_pad(left, &left_padded)) {
                SIFT3D_ERR("draw_matches: unable to pad images \n");
                return SIFT3D_FAILURE;
	}

	// Draw a concatenated image
	if (im_concat(&left_padded, &right_padded, 0, concat_arg)) {
                SIFT3D_ERR("draw_matches: Could not concatenate the "
                        "images \n");
                goto draw_matches_quit;
        }

        // Optionally draw the keypoints
        if (keys != NULL) { 

                // Convert inputs to double
                if (convert_Mat_rm(keys_right, &keys_right_draw, 
			SIFT3D_DOUBLE) ||
                        convert_Mat_rm(keys_left, &keys_left_draw, 
				SIFT3D_DOUBLE))
                        goto draw_matches_quit;
       
                // Pad the x-coordinate 
                for (i = 0; i < keys_right->num_rows; i++) {
                        SIFT3D_MAT_RM_GET(&keys_right_draw, i, 0, double) += 
                                right_pad; 
                }

                // Concatenate the points
                if (concat_Mat_rm(&keys_left_draw, &keys_right_draw,
                        &keys_draw, 0))
                        goto draw_matches_quit;

                // Draw the points
                if (draw_points(&keys_draw, SIFT3D_IM_GET_DIMS(concat_arg), 
                        1, keys))
                        goto draw_matches_quit;
        }

	// Optionally draw the lines 
        if (lines != NULL) {

                // Convert input to double
                if (convert_Mat_rm(match_right, &match_right_draw, 
			SIFT3D_DOUBLE))
                        goto draw_matches_quit;

                // Pad the x-coordinate 
                for (i = 0; i < match_right->num_rows; i++) {
                        SIFT3D_MAT_RM_GET(&match_right_draw, i, 0, double) += 
                                right_pad;
                }

                // Draw the lines
                if (draw_lines(match_left, &match_right_draw, 
                        SIFT3D_IM_GET_DIMS(concat_arg), lines))
                        goto draw_matches_quit;
        }

        // Clean up
        im_free(&concat_temp);
        im_free(&left_padded);
        im_free(&right_padded); 
        cleanup_Mat_rm(&keys_right_draw); 
        cleanup_Mat_rm(&keys_left_draw); 
        cleanup_Mat_rm(&keys_draw); 
        cleanup_Mat_rm(&match_right_draw); 
	return SIFT3D_SUCCESS;

draw_matches_quit:
        im_free(&concat_temp);
        im_free(&left_padded);
        im_free(&right_padded); 
        cleanup_Mat_rm(&keys_right_draw);
        cleanup_Mat_rm(&keys_left_draw); 
        cleanup_Mat_rm(&keys_draw); 
        cleanup_Mat_rm(&match_right_draw);
        return SIFT3D_FAILURE;
}

/* Write a Keypoint_store to a text file. The keypoints are stored in a matrix
 * (.csv, .csv.gz), where each keypoint is a row. The elements of each row are
 * as follows:
 *
 * x y z o s ori11 ori12 ... orinn
 *
 * x - the x-coordinate
 * y - the y-coordinate
 * z - the z-coordinate
 * o - the pyramid octave. To convert to image coordinates, multiply x,y,z by 
 *      pow(2, o)
 * s - the scale coordinate
 * ori(ij) - the ith row, jth column of the orientation matrix */
int write_Keypoint_store(const char *path, const Keypoint_store *const kp) {

        Mat_rm mat;
	int i, i_R, j_R;

        // Keypoint data format constants
        const int kp_x = 0; // column of x-coordinate
        const int kp_y = 1; // column of y-coordinate
        const int kp_z = 2; // column of z-coordinate
        const int kp_o = 3; // column of octave index
        const int kp_s = 4; // column of s-coordinate
        const int kp_ori = 5; // first column of the orientation matrix
        const int ori_numel = IM_NDIMS * IM_NDIMS; // Number of orientation 
                // elements
        const int num_rows = kp->slab.num;
        const int num_cols = kp_ori + ori_numel;

        // Initialize the matrix
        if (init_Mat_rm(&mat, num_rows, num_cols, SIFT3D_DOUBLE, 
		SIFT3D_FALSE))
                return SIFT3D_FAILURE;
       
        // Write the keypoints 
        for (i = 0; i < num_rows; i++) {

                const Keypoint *const key = kp->buf + i;
                const Mat_rm *const R = &key->R;

                // Write the coordinates 
                SIFT3D_MAT_RM_GET(&mat, i, kp_x, double) = key->xd;
                SIFT3D_MAT_RM_GET(&mat, i, kp_y, double) = key->yd; 
                SIFT3D_MAT_RM_GET(&mat, i, kp_z, double) = key->zd; 
                SIFT3D_MAT_RM_GET(&mat, i, kp_o, double) = key->o; 
                SIFT3D_MAT_RM_GET(&mat, i, kp_s, double) = key->sd;

                // Write the orientation matrix
                SIFT3D_MAT_RM_LOOP_START(R, i_R, j_R)

                        const int kp_idx = kp_ori + 
                                SIFT3D_MAT_RM_GET_IDX(R, i_R, j_R);
        
                        SIFT3D_MAT_RM_GET(&mat, i, kp_idx, double) = 
                                (double) SIFT3D_MAT_RM_GET(R, i_R, j_R, float);

                SIFT3D_MAT_RM_LOOP_END
        }

        // Write the matrix 
        if (write_Mat_rm(path, &mat))
                goto write_kp_quit;

        // Clean up
        cleanup_Mat_rm(&mat);

        return SIFT3D_SUCCESS;

write_kp_quit:
        cleanup_Mat_rm(&mat);
        return SIFT3D_FAILURE;
}

/* Write SIFT3D descriptors to a text file.
 * See SIFT3D_Descriptor_store_to_Mat_rm for the file format. */
int write_SIFT3D_Descriptor_store(const char *path, 
        const SIFT3D_Descriptor_store *const desc) {

        Mat_rm mat;

        // Initialize the matrix
        if (init_Mat_rm(&mat, 0, 0, SIFT3D_FLOAT, SIFT3D_FALSE))
                return SIFT3D_FAILURE;
     
        // Write the data into the matrix 
        if (SIFT3D_Descriptor_store_to_Mat_rm(desc, &mat))
                goto write_desc_quit;

        // Write the matrix to the file
        if (write_Mat_rm(path, &mat))
                goto write_desc_quit;

        // Clean up
        cleanup_Mat_rm(&mat);
        return SIFT3D_SUCCESS;

write_desc_quit:
        cleanup_Mat_rm(&mat);
        return SIFT3D_FAILURE;
}

