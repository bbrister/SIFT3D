/* sift.c
 * ----------------------------------------------------------------
 * Rice MRI Team
 * ----------------------------------------------------------------
 * This file contains all routines needed to initialize, delete, 
 * and run the 3D SIFT detector and descriptor.
 *-----------------------------------------------------------------
 * Created: Blaine Rister 12/26/2013
 * Last updated: Blaine Rister 12/13/2014
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <getopt.h>
#include "types.h"
#include "macros.h"
#include "imutil.h"

#include "sift.h"

/* Internal parameters */
#define ORI_LB 0.1		// Lower bound ratio on orientation magnitude
#define ORI_UB 0.95		// Upper bound ratio on orientation magnitude
#define ORI_GRAD_THRESH 1E-10   // Minimum norm of average gradient
//#define EIG_MIN 1E-2		// Minimum allowed eigenvalue magnitude
#define EIG_MIN 1E-10		// Minimum allowed eigenvalue magnitude
#define EIG_MAX_RATIO 0.90	// Maximum ratio of eigenvalue magnitudes
#define EIG_COS_ANGLE_MIN 0.5	// Minimum |cos(angle)| between eigenvector and gradient
#define BARY_EPS FLT_EPSILON * 1E1	// Error tolerance for barycentric coordinates
#define ORI_SIG_FCTR 1.5        // Ratio of window parameter to keypoint scale
#define ORI_RAD_FCTR 3.0 // Ratio of window radius to parameter
#define DESC_SIG_FCTR 7.071067812 // See, ORI_SIG_FCTR, 5 * sqrt(2)
#define DESC_RAD_FCTR 2.0  // See ORI_RAD_FCTR
#define TRUNC_THRESH (0.2f * 128.0f / DESC_NUMEL) // Descriptor truncation threshold

/* Internal return codes */
#define REJECT 1

/* Internal constants */
#define GR 1.6180339887	// Golden ratio
#define PRECOMP_NC 4    // Number of channels in precomputed gradient image

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
 * Delimit with IM_LOOP_END. */
#define IM_LOOP_SPHERE_START(im, x, y, z, vcenter, rad, vdisp, sq_dist) \
	const int x_start = (int) MAX((int) (vcenter)->x - (int) ((rad) + 0.5), 1); \
	const int x_end   = (int) MIN((int) (vcenter)->x + (int) ((rad) + 0.5), im->nx - 2); \
	const int y_start = (int) MAX((int) (vcenter)->y - (int) ((rad) + 0.5), 1); \
	const int y_end   = (int) MIN((int) (vcenter)->y + (int) ((rad) + 0.5), im->ny - 2); \
	const int z_start = (int) MAX((int) (vcenter)->z - (int) ((rad) + 0.5), 1); \
	const int z_end   = (int) MIN((int) (vcenter)->z + (int) ((rad) + 0.5), im->nz - 2); \
	IM_LOOP_LIMITED_START(im, x, y, z, x_start, x_end, y_start, y_end, \
			      z_start, z_end) \
	    (vdisp)->x = ((float) x + 0.5f) - (vcenter)->x; \
	    (vdisp)->y = ((float) y + 0.5f) - (vcenter)->y; \
	    (vdisp)->z = ((float) z + 0.5f) - (vcenter)->z; \
	    (sq_dist) = CVEC_L2_NORM_SQ(vdisp); \
	    if ((sq_dist) > (rad) * (rad)) \
		continue; \

/* Global variables */
extern CL_data cl_data;

/* Helper routines */
static int init_geometry(SIFT3D *sift3d);
static int resize_SIFT3D(SIFT3D *sift3d);
static int build_gpyr(SIFT3D *sift3d);
static int build_dog(SIFT3D *dog);
static int detect_extrema(SIFT3D *sift3d, Keypoint_store *kp);
static int refine_keypoints(SIFT3D *sift3d, Keypoint_store *kp);
static int assign_hist_ori(const Image *const im, const Cvec *const vcenter,
                          const double sigma, Mat_rm *const R);
static int assign_eig_ori(const Image *const im, const Cvec *const vcenter,
                          const double sigma, Mat_rm *const R);
static int assign_orientations(SIFT3D *sift3d, Keypoint_store *kp);
static int Cvec_to_sbins(const Cvec * const vd, Svec * const bins);
static void refine_Hist(Hist *hist);
static int init_cl_SIFT3D(SIFT3D *sift3d);
static int cart2bary(const Cvec * const cart, const Tri * const tri, 
		      Cvec * const bary, float * const k);
static void SIFT3D_desc_acc_interp(const SIFT3D * const sift3d, 
				   const Cvec * const vbins, 
				   const Cvec * const grad,
				   SIFT3D_Descriptor * const desc);
static void extract_descrip(SIFT3D *sift3d, Image *im,
	   Keypoint *key, SIFT3D_Descriptor *desc);
static int argv_permute(const int argc, char *const *argv, 
                        const unsigned char *processed, 
                        const int optind_start);
static int extract_dense_descriptors(SIFT3D *const sift3d,
        const Image *const in, Image *const desc);
static int extract_dense_descriptors_rotate(SIFT3D *const sift3d,
        const Image *const in, Image *const desc);
static int vox2hist(const Image *const im, const int x, const int y,
        const int z, Hist *const hist);
static int hist2vox(Hist *const hist, const Image *const im, const int x, 
        const int y, const int z);

/* Initialize geometry tables. */
static int init_geometry(SIFT3D *sift3d) {

	Mat_rm V, F;
	Cvec temp1, temp2, temp3, n;
	float mag;
	int i, j;

	Mesh * const mesh = &sift3d->mesh;

	/* Verices of a regular icosahedron inscribed in the unit sphere. */
	const float vert[] = {  0,  1,  GR,
			        0, -1,  GR,
			        0,  1, -GR,
			        0, -1, -GR,
			        1,  GR,  0,
			       -1,  GR,  0,
			        1, -GR,  0,
			       -1, -GR,  0,
			       GR,   0,  1,
			      -GR,   0,  1,
			       GR,   0, -1, 
			      -GR,   0, -1 }; 

	/* Vertices of the faces of the icosahedron. */
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
	if (init_Mat_rm_p(&V, vert, ICOS_NVERT, 3, FLOAT, FALSE) ||
	    init_Mat_rm_p(&F, faces, ICOS_NFACES, 3, FLOAT, FALSE))
		return FAILURE;
			    
	// Initialize triangle memory
	if ((mesh->tri = (Tri *) malloc(ICOS_NFACES * sizeof(Tri))) == NULL)
		return FAILURE;
 
	// Populate the triangle struct for each face
	for (i = 0; i < ICOS_NFACES; i++) {

		Tri * const tri = mesh->tri + i;	
		Cvec * const v = tri->v;

		// Initialize the vertices
		for (j = 0; j < 3; j++) {

			const float mag_expected = sqrt(1 + GR * GR);
			int * const idx = tri->idx + j;

			*idx = MAT_RM_GET(&F, i, j, float);

			// Initialize the vector
			v[j].x = MAT_RM_GET(&V, *idx, 0, float);
			v[j].y = MAT_RM_GET(&V, *idx, 1, float);
			v[j].z = MAT_RM_GET(&V, *idx, 2, float);

			// Normalize to unit length
			mag = CVEC_L2_NORM(v + j);
			assert(fabsf(mag - mag_expected) < 1E-10);
			CVEC_SCALE(v + j, 1.0f / mag);
		}

		// Compute the normal vector at v[0] as  (V2 - V1) X (V1 - V0)
		CVEC_OP(v + 2, v + 1, -, &temp1);
		CVEC_OP(v + 1, v, -, &temp2);
		CVEC_CROSS(&temp1, &temp2, &n);

		// Ensure this vector is facing outward from the origin
		if (CVEC_DOT(&n, v) < 0) {
			// Swap two vertices
			temp1 = v[0];
			v[0] = v[1];
			v[1] = temp1;

			// Compute the normal again
			CVEC_OP(v + 2, v + 1, -, &temp1);
			CVEC_OP(v + 1, v, -, &temp2);
			CVEC_CROSS(&temp1, &temp2, &n);
		}
		assert(CVEC_DOT(&n, v) >= 0);

		// Ensure the triangle is equilateral
		CVEC_OP(v + 2, v, -, &temp3);
		assert(fabsf(CVEC_L2_NORM(&temp1) - CVEC_L2_NORM(&temp2)) < 
			1E-10);
		assert(fabsf(CVEC_L2_NORM(&temp1) - CVEC_L2_NORM(&temp3)) < 
			1E-10);
	}	
	
	return SUCCESS;
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

	CVEC_OP(v + 1, v, -, &e1);
	CVEC_OP(v + 2, v, -, &e2);
	CVEC_CROSS(cart, &e2, &p);
	det = CVEC_DOT(&e1, &p);

	// Reject unstable points
	if (fabsf(det) < BARY_EPS) {
		return FAILURE;
	}

	det_inv = 1.0f / det;

	t = v[0];
	CVEC_SCALE(&t, -1.0f);	

	CVEC_CROSS(&t, &e1, &q);

	bary->y = det_inv * CVEC_DOT(&t, &p);	
	bary->z = det_inv * CVEC_DOT(cart, &q);
	bary->x = 1.0f - bary->y - bary->z;

	*k = CVEC_DOT(&e2, &q) * det_inv;

#ifndef NDEBUG
	Cvec temp1, temp2, temp3;
        double residual;

	// Verify k * c = bary->x * v1 + bary->y * v2 + bary->z * v3
	temp1 = v[0];
	temp2 = v[1];
	temp3 = v[2];
	CVEC_SCALE(&temp1, bary->x);
	CVEC_SCALE(&temp2, bary->y);	
	CVEC_SCALE(&temp3, bary->z);	
	CVEC_OP(&temp1, &temp2, +, &temp1);
	CVEC_OP(&temp1, &temp3, +, &temp1);
	CVEC_SCALE(&temp1, 1.0f / *k);
	CVEC_OP(&temp1, cart, -, &temp1);
        residual = CVEC_L2_NORM(&temp1);
	if (residual > BARY_EPS) {
                printf("cart2bary: residual: %f", residual);
                exit(1);
        }
#endif
	return SUCCESS;
}

/* Briefly initialize a Keypoint_store for first use.
 * This does not need to be called to re-use the store
 * for a new image. */
void init_Keypoint_store(Keypoint_store *kp) {
	init_Slab(&kp->slab);
	kp->buf = (Keypoint *) kp->slab.buf;
}

/* Briefly initialize a SIFT_Descriptor_store for first use.
* This does not need to be called to re-use the store
* for a new image. */
void init_SIFT3D_Descriptor_store(SIFT3D_Descriptor_store *desc) {
	desc->buf = NULL;
}

/* Initializes the OpenCL data for this SIFT3D struct. This
 * increments the reference counts for shared data. */
static int init_cl_SIFT3D(SIFT3D *sift3d) {
#ifdef USE_OPENCL
	cl_image_format image_format;

	// Initialize basic OpenCL platform and context info
	image_format.image_channel_order = CL_R;
	image_format.image_channel_data_type = CL_FLOAT;
	if (init_cl(&cl_data, PLATFORM_NAME_NVIDIA, CL_DEVICE_TYPE_GPU,
 		    CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, 
                    image_format))
		return FAILURE;

	// Load and compile the downsampling kernel

#endif
	return SUCCESS;
}

void set_first_octave_SIFT3D(SIFT3D *const sift3d, 
                                const int first_octave) {
	sift3d->dog.first_octave = sift3d->gpyr.first_octave = first_octave;
}

void set_peak_thresh_SIFT3D(SIFT3D *const sift3d,
                                const double peak_thresh) {
        sift3d->peak_thresh = peak_thresh;
}

void set_corner_thresh_SIFT3D(SIFT3D *const sift3d,
                                const double corner_thresh) {
        sift3d->corner_thresh = corner_thresh;
}

void set_num_octaves_SIFT3D(SIFT3D *const sift3d,
                                const int num_octaves) {
	sift3d->dog.num_octaves = sift3d->gpyr.num_octaves = num_octaves;
}

void set_num_kp_levels_SIFT3D(SIFT3D *const sift3d,
                                const int num_kp_levels) {

	const int num_dog_levels = num_kp_levels + 2;
	const int num_gpyr_levels = num_dog_levels + 1;

	sift3d->dog.num_kp_levels = sift3d->gpyr.num_kp_levels = 
                num_kp_levels;
	sift3d->dog.num_levels = num_dog_levels;
	sift3d->gpyr.num_levels = num_gpyr_levels;
}

void set_sigma_n_SIFT3D(SIFT3D *const sift3d,
                                const double sigma_n) {
	sift3d->dog.sigma_n = sift3d->gpyr.sigma_n = sigma_n;
}

void set_sigma0_SIFT3D(SIFT3D *const sift3d,
                                const double sigma0) {
	sift3d->dog.sigma0 = sift3d->gpyr.sigma0 = sigma0;
}

/* Initialize a SIFT3D struct with the default parameters. */
int init_SIFT3D(SIFT3D *sift3d) {

	int num_dog_levels, num_gpyr_levels;

	// Initialize to defaults
	const int first_octave = DEFAULT_FIRST_OCTAVE;
	const double peak_thresh = DEFAULT_PEAK_THRESH;
	const double corner_thresh = DEFAULT_CORNER_THRESH;
	const int num_octaves = -1;
	const int num_kp_levels = DEFAULT_NUM_KP_LEVELS;
	const double sigma_n = DEFAULT_SIGMA_N;
	const double sigma0 = DEFAULT_SIGMA0;
        const int dense_rotate = FALSE;

	// First-time pyramid initialization
	sift3d->dog.levels = sift3d->gpyr.levels = NULL;

        // Intialize the geometry tables
	if (init_geometry(sift3d))
		return FAILURE;

	// init static OpenCL programs and contexts, if support is enabled
	if (init_cl_SIFT3D(sift3d))
		return FAILURE;

	// Initialize image to null, to mark for resizing
	sift3d->im = NULL;

	// Declare pyramid memory null, to mark for initialization 
        // (see init_pyramid)
	sift3d->gpyr.levels = sift3d->dog.levels = NULL;

	// Save data
	sift3d->dog.first_level = sift3d->gpyr.first_level = -1;
        set_first_octave_SIFT3D(sift3d, first_octave);
        set_num_octaves_SIFT3D(sift3d, num_octaves);
        set_num_kp_levels_SIFT3D(sift3d, num_kp_levels);
        set_sigma_n_SIFT3D(sift3d, sigma_n);
        set_sigma0_SIFT3D(sift3d, sigma0);
        set_peak_thresh_SIFT3D(sift3d, peak_thresh);
        set_corner_thresh_SIFT3D(sift3d, corner_thresh);
        sift3d->dense_rotate = dense_rotate;

	return SUCCESS;
}

/* Helper function to permute the arguments so that all unprocessed arguments
 * are at the end of the array. This function only processes the closed 
 * interval [optind_start, infinity]. The variable optind is reset to the
 * first unprocessed argument in this interval, in the final argv array. */
static int argv_permute(const int argc, char *const *argv, 
                        const unsigned char *processed, 
                        const int optind_start) {

        char **next_processed, **next_unprocessed;
        char **argv_in;
        int i, num_processed;

        const size_t argv_size = argc * sizeof(char *);

        // Allocate memory
        if ((argv_in = (char **) malloc(argv_size)) == NULL) {
                fprintf(stderr, "argv_permute: out of memory \n");
                return FAILURE;
        }

        // Copy the input argv array 
        memcpy(argv_in, argv, argv_size);

        // Get the number of processed arguments
        num_processed = 0;
        for (i = optind_start; i < argc; i++) {
                if (!processed[i])
                        continue;
                num_processed++;
        }

        // Get the starting locations for each section of argv
        next_processed = (char **) argv + optind_start;
        next_unprocessed = next_processed + num_processed;

        // Build the permuted argv in the place of the original 
        for (i = optind_start; i < argc; i++) {
                
                char *const arg = argv_in[i];

                if (processed[i]) {
                        *next_processed++ = arg;
                } else {
                        *next_unprocessed++ = arg;
                } 
        }

        // Clean up
        free(argv_in);

        // Reset optind
        optind = optind_start + num_processed;

        return SUCCESS;
}

/* Set the parameters of a SIFT3D struct from the given command line 
 * arguments. The argument SIFT3D must be initialized with
 * init_SIFT3D prior to calling this function. All options not
 * specified in argv will remain at their previous values.
 *
 * NOTE: This function uses the POSIX convention for getopt, i.e. processing 
 * stops at the first non-option argument. Due to a bug in some getopt
 * implementations, any subsequent calls to getopt or getopt_long will
 * follow the POSIX convention, regardless of the arguments you provide. 
 * Thus, any process calling this funciton must adopt the POSIX convention.
 *
 * Options:
 * --first_octave	 - the first octave (int)
 * --peak_thresh	 - threshold on DoG extrema magnitude (double)
 * --corner_thresh - threshold on edge score (double)
 * --num_octaves	 - total number of octaves (default: process
 *	  		until one dimension is less than 8, use 
 *			-1 for default) (int)
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
 *      optind_ret - If NULL, does nothing. Otherwise, an index to the remaining
 *             options is written to this address, as in GNU getopt.
 *      check_err - If nonzero, report unrecognized options as errors
 * Return value: 
 *      0 on success, -1 on error. */
int parse_args_SIFT3D(SIFT3D *const sift3d,
        const int argc, char *const *argv, int *optind_ret, 
        const int check_err) {

        unsigned char *processed;
        double dval;
        int c, err, ival;

        // Options
        const struct option longopts[] = {
                {"first_octave", required_argument, NULL, 'a'},
                {"peak_thresh", required_argument, NULL, 'b'},
                {"corner_thresh", required_argument, NULL, 'c'},
                {"num_octaves", required_argument, NULL, 'd'},
                {"num_kp_levels", required_argument, NULL, 'e'},
                {"sigma_n", required_argument, NULL, 'f'},
                {"sigma0", required_argument, NULL, 'g'}
        };

        // Starting getopt variables 
        const int optind_start = optind;
        const int opterr_start = opterr;

        // Set the error checking behavior
        opterr = check_err;

        // Intialize intermediate data
        if ((processed = calloc(argc, sizeof(char *))) == NULL) {
                fprintf(stderr, "parse_args_SIFT3D: out of memory\n");
                return FAILURE;
        }
        err = FALSE;

        // Process the arguments
        while ((c = getopt_long(argc, argv, "+", longopts, NULL)) != -1) {

                const int idx = optind - 1;

                // Convert the value to double and integer
                if (optarg != NULL) {
                        dval = atof(optarg);
                        ival = atoi(optarg);
                }

                switch (c) {
                        case 'a':
                                set_first_octave_SIFT3D(sift3d, 
                                        ival);
                                processed[idx - 1] = TRUE;
                                processed[idx] = TRUE;
                                break;
                        case 'b':
                                set_peak_thresh_SIFT3D(sift3d, 
                                        dval);
                                processed[idx - 1] = TRUE;
                                processed[idx] = TRUE;
                                break;
                        case 'c':
                                set_corner_thresh_SIFT3D(sift3d, 
                                        dval);
                                processed[idx - 1] = TRUE;
                                processed[idx] = TRUE;
                                break;
                        case 'd':
                                set_num_octaves_SIFT3D(sift3d, 
                                        ival);
                                processed[idx - 1] = TRUE;
                                processed[idx] = TRUE;
                                break;
                        case 'e':
                                set_num_kp_levels_SIFT3D(sift3d, 
                                        ival);
                                processed[idx - 1] = TRUE;
                                processed[idx] = TRUE;
                                break;
                        case 'f':
                                set_sigma_n_SIFT3D(sift3d, dval);
                                processed[idx - 1] = TRUE;
                                processed[idx] = TRUE;
                                break;
                        case 'g':
                                set_sigma0_SIFT3D(sift3d, dval);
                                processed[idx - 1] = TRUE;
                                processed[idx] = TRUE;
                                break;
                        case '?':
                        default:
                                if (!check_err)
                                        break;
                                err = TRUE;
                }
        }

        // Put all unprocessed options at the end
        if (argv_permute(argc, argv, processed, optind_start))
                goto parse_args_quit;

        // Return to the default settings
        opterr = opterr_start;

        // Return optind
        if (optind_ret != NULL)
                *optind_ret = optind;

        // Clean up
        free(processed);

        // Return the error condition, if error checking is enabled
        if (check_err && err)
                return FAILURE;

        return SUCCESS;

parse_args_quit:
        free(processed);
        return FAILURE;
}

/* Helper routine to resize a SIFT3D struct and recompile filters, 
 * according to the dimensions of im. */
static int resize_SIFT3D(SIFT3D *sift3d) {

	double nx, ny, nz;
	int last_octave, num_octaves; 

	Image const *im = sift3d->im;
	const int first_octave = sift3d->gpyr.first_octave;

	// Compute the number of octaves, if not specified by user
	if ((num_octaves = sift3d->gpyr.num_octaves) == -1) {
		last_octave = (int) log2((double) MIN(MIN(im->nx, im->ny),
			im->nz)) - 3 - first_octave;		// min size: 8 in any dimension

		num_octaves = last_octave - first_octave + 1;
	} else {
		// Number of octaves specified by user: compute last octave
		last_octave = num_octaves + first_octave - 1;
	}

	// Update pyramid data
	sift3d->gpyr.last_octave = sift3d->dog.last_octave = last_octave;
	sift3d->dog.num_octaves = sift3d->gpyr.num_octaves = num_octaves;

	// Re-initialize the pyramid
	if (init_pyramid(&sift3d->gpyr, sift3d->im) ||
		init_pyramid(&sift3d->dog, sift3d->im))
		return FAILURE;

	// Compute the Gaussian filters
	if (init_gss(&sift3d->filters.gss, &sift3d->gpyr))
		return FAILURE;

	return SUCCESS;
}

/* Build the GSS pyramid on a single CPU thread */
static int build_gpyr(SIFT3D *sift3d) {

	Sep_FIR_filter *f;
	Image *cur, *prev;
	int o, s;

	Pyramid *const gpyr = &sift3d->gpyr;
	const GSS_filters *const gss = &sift3d->filters.gss;
	const int s_start = gpyr->first_level + 1;
	const int s_end = gpyr->last_level;
	const int o_start = gpyr->first_octave;
	const int o_end = gpyr->last_octave;

	// Build the first image
	cur = PYR_IM_GET(gpyr, o_start, s_start - 1);
	prev = sift3d->im;
#ifdef USE_OPENCL
	if (im_load_cl(cur, FALSE))
		return FAILURE;	
#endif

	f = (Sep_FIR_filter *) &gss->first_gauss.f;
	if (apply_Sep_FIR_filter(prev, cur, f))
		return FAILURE;

	// Build the rest of the pyramid
	PYR_LOOP_LIMITED_START(o, s, o_start, o_end, s_start, s_end)
			cur = PYR_IM_GET(gpyr, o, s);
			prev = PYR_IM_GET(gpyr, o, s - 1);
			f = &gss->gauss_octave[s].f;
			if (apply_Sep_FIR_filter(prev, cur, f))
				return FAILURE;
#ifdef USE_OPENCL
			if (im_read_back(cur, FALSE))
				return FAILURE;
#endif
		PYR_LOOP_SCALE_END
		// Downsample
		if (o != o_end) {
			prev = cur;
			cur = PYR_IM_GET(gpyr, o + 1, s_start - 1);
			if (im_downsample_2x(prev, cur))
				return FAILURE;
		}
	PYR_LOOP_OCTAVE_END

#ifdef USE_OPENCL
	clFinish_all();
#endif

	return SUCCESS;
}

static int build_dog(SIFT3D *sift3d) {

	Image *gpyr_cur, *gpyr_next, *dog_level;
	int o, s;

	Pyramid *const dog = &sift3d->dog;
	Pyramid *const gpyr = &sift3d->gpyr;

	PYR_LOOP_START(dog, o, s)
		gpyr_cur = PYR_IM_GET(gpyr, o, s);
		gpyr_next = PYR_IM_GET(gpyr, o, s + 1);			
		dog_level = PYR_IM_GET(dog, o, s);
		
		if (im_subtract(gpyr_cur, gpyr_next, 
						dog_level))
			return FAILURE;
	PYR_LOOP_END

	return SUCCESS;
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
	const int o_end = dog->last_octave;
	const int s_start = dog->first_level + 1;
	const int s_end = dog->last_level - 1;

	// Verify the inputs
	if (dog->num_levels < 3) {
		printf("detect_extrema: Requires at least 3 levels per octave, "
			   "proivded only %d", dog->num_levels);
		return FAILURE;
	}

	// Initialize dimensions of keypoint store
	cur = PYR_IM_GET(dog, o_start, s_start);
	kp->nx = cur->nx;
	kp->ny = cur->ny;
	kp->nz = cur->nz;

#define CMP_NEIGHBORS(im, x, y, z, CMP, IGNORESELF, val) ( \
	(val) CMP IM_GET_VOX( (im), (x),     (y),    (z) - 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) - 1, (y),    (z) - 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) + 1, (y),    (z) - 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x),     (y) - 1,(z) - 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x),     (y) + 1,(z) - 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) - 1, (y) - 1,(z) - 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) + 1, (y) - 1,(z) - 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) - 1, (y) + 1,(z) - 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) + 1, (y) + 1,(z) - 1, 0) && \
	((val) CMP IM_GET_VOX( (im), (x),     (y),    (z), 0   ) || \
	    IGNORESELF) && \
	(val) CMP IM_GET_VOX( (im), (x) - 1, (y),    (z), 0    ) && \
	(val) CMP IM_GET_VOX( (im), (x) + 1, (y),    (z), 0    ) && \
	(val) CMP IM_GET_VOX( (im), (x),     (y) - 1,(z), 0    ) && \
	(val) CMP IM_GET_VOX( (im), (x),     (y) + 1,(z), 0    ) && \
	(val) CMP IM_GET_VOX( (im), (x) - 1, (y) - 1,(z), 0    ) && \
	(val) CMP IM_GET_VOX( (im), (x) + 1, (y) - 1,(z), 0    ) && \
	(val) CMP IM_GET_VOX( (im), (x) - 1, (y) + 1,(z), 0    ) && \
	(val) CMP IM_GET_VOX( (im), (x) + 1, (y) + 1,(z), 0    ) && \
	(val) CMP IM_GET_VOX( (im), (x),     (y),    (z) + 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) - 1, (y),    (z) + 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) + 1, (y),    (z) + 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x),     (y) - 1,(z) + 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x),     (y) + 1,(z) + 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) - 1, (y) - 1,(z) + 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) + 1, (y) - 1,(z) + 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) - 1, (y) + 1,(z) + 1, 0) && \
	(val) CMP IM_GET_VOX( (im), (x) + 1, (y) + 1,(z) + 1, 0) )

	num = 0;
	PYR_LOOP_LIMITED_START(o, s, o_start, o_end, s_start, s_end)  

		// Select current and neighboring levels
		prev = PYR_IM_GET(dog, o, s - 1);
		cur = PYR_IM_GET(dog, o, s);
		next = PYR_IM_GET(dog, o, s + 1);

		// Find maximum DoG value at this level
		dogmax = 0.0f;
		IM_LOOP_START(cur, x, y, z)
			dogmax = MAX(dogmax, 
                                fabsf(IM_GET_VOX(cur, x, y, z, 0)));
		IM_LOOP_END
		// Adjust threshold
		peak_thresh = sift3d->peak_thresh * dogmax;

		// Loop through all non-boundary pixels
		x_start = y_start = z_start = 1;
		x_end = cur->nx - 2;
		y_end = cur->ny - 2;
		z_end = cur->nz - 2;
		IM_LOOP_LIMITED_START(cur, x, y, z, x_start, x_end, y_start,
							  y_end, z_start, z_end)
			// Sample the center value
			pcur = IM_GET_VOX(cur, x, y, z, 0);

			// Apply the peak threshold
			if ((pcur > peak_thresh || pcur < -peak_thresh) && ((
				// Compare to the neighbors
				CMP_NEIGHBORS(prev, x, y, z, >, FALSE, pcur) &&
				CMP_NEIGHBORS(cur, x, y, z, >, TRUE, pcur) &&
				CMP_NEIGHBORS(next, x, y, z, >, FALSE, pcur)
				) || (
				CMP_NEIGHBORS(prev, x, y, z, <, FALSE, pcur) &&
				CMP_NEIGHBORS(cur, x, y, z, <, TRUE, pcur) &&
				CMP_NEIGHBORS(next, x, y, z, <, FALSE, pcur))))
				{
					// Add a keypoint candidate
					num++;
					RESIZE_KP_STORE(kp, num, sizeof(Keypoint));
					key = kp->buf + num - 1;
					key->o = o;
					key->s = s;
					key->xi = x;
					key->yi = y;
					key->zi = z;

				}
		IM_LOOP_END
	PYR_LOOP_END
#undef CMP_NEIGHBORS

	return SUCCESS;
}

/* Refine keypoint locations to sub-pixel accuracy. */
static int refine_keypoints(SIFT3D *sift3d, Keypoint_store *kp) {

	Mat_rm B, Hi, Hs, X;
IGNORE_UNUSED
	Cvec vd;
	double xd, yd, zd, sd;
IGNORE_UNUSED
	int x, y, z, c, xnew, ynew, znew, i, j, k, l;

	// Initialize data 
	init_Mat_rm(&B, 4, 1, DOUBLE, FALSE);
	init_Mat_rm(&X, 4, 1, DOUBLE, FALSE);
	init_Mat_rm(&Hi, 3, 3, DOUBLE, FALSE);
	init_Mat_rm(&Hs, 4, 4, DOUBLE, FALSE);

	for (k = 0; k < kp->slab.num; k++) {

		// Misc. constant data
		Keypoint * const key = kp->buf + k;
		const int o = key->o;
		const int s = key->s;
		const Image const *prev = PYR_IM_GET(&sift3d->dog, o, s - 1);
		const Image const *cur = PYR_IM_GET(&sift3d->dog, o, s);
		const Image const *next = PYR_IM_GET(&sift3d->dog, o, s + 1);

		// Bound the translation to all non-boundary pixels
		const double xmin = 1;
		const double ymin = 1;
		const double zmin = 1;
		const double xmax = cur->nx - 2 - DBL_EPSILON;
		const double ymax = cur->ny - 2 - DBL_EPSILON;
		const double zmax = cur->nz - 2 - DBL_EPSILON;
	    
		// Bound the scale to that of the neighboring levels
		const double smin = prev->s;
		const double smax = next->s;
	
		// Initialize mutable data	
		x = key->xi;
		y = key->yi;
		z = key->zi;
		xd = (double) x + 0.5;
		yd = (double) y + 0.5;
		zd = (double) z + 0.5;
		sd = cur->s; 

		// Refine the keypoint for a fixed number of iterations
		for (l = 0; l < 5; l++) {

		assert(x >= 1 && y >= 1 && z >= 1 && x <= cur->nx - 2 &&
		       y <= cur->ny -2 && z <= cur->nz - 2); 

#define PARABOLA
#ifndef PARABOLA 
		// Form the gradient
		IM_GET_GRAD(cur, x, y, z, 0, &vd);

		// Form the response vector as the negative gradient
		MAT_RM_GET(&B, 0, 0, double) = (double) -vd.x;
		MAT_RM_GET(&B, 1, 0, double) = (double) -vd.y;
		MAT_RM_GET(&B, 2, 0, double) = (double) -vd.z;
		MAT_RM_GET(&B, 3, 0, double) = 
		   (double) - 0.5 * (IM_GET_VOX(next, x, y, z, 0) - 
			      IM_GET_VOX(prev, x, y, z, 0));

		// Form the Hessian
		IM_GET_HESSIAN(cur, x, y, z, c, 0, &Hi, double);
		MAT_RM_LOOP_START(&Hi, i, j)
			MAT_RM_GET(&Hs, i, j, double) = 
				MAT_RM_GET(&Hi, i, j, double);
		MAT_RM_LOOP_END

		// Dsx
		MAT_RM_GET(&Hs, 0, 3, double) = 
		MAT_RM_GET(&Hs, 3, 0, double) =
			(double) 0.25 * (IM_GET_VOX(next, x + 1, y, z, 0) -
			 IM_GET_VOX(prev, x + 1, y, z, 0) + 
			 IM_GET_VOX(prev, x - 1, y, z, 0) - 
			 IM_GET_VOX(next, x - 1, y, z, 0)); 

		// Dsy 
		MAT_RM_GET(&Hs, 1, 3, double) = 
		MAT_RM_GET(&Hs, 3, 1, double) = 
			(double) 0.25 * (IM_GET_VOX(next, x, y + 1, z, 0) -
			 IM_GET_VOX(prev, x, y + 1, z, 0) + 
			 IM_GET_VOX(prev, x, y - 1, z, 0) - 
			 IM_GET_VOX(next, x, y - 1, z, 0)); 

		// Dsz 
		MAT_RM_GET(&Hs, 2, 3, double) = 
		MAT_RM_GET(&Hs, 3, 2, double) = 
			(double) 0.25 * (IM_GET_VOX(next, x, y, z + 1) -
			 IM_GET_VOX(prev, x, y, z + 1) + 
			 IM_GET_VOX(prev, x, y, z - 1) - 
			 IM_GET_VOX(next, x, y, z - 1)); 

		// Dss  
		MAT_RM_GET(&Hs, 3, 3, double) = 
			(double) 0.25 * (IM_GET_VOX(next, x, y, z, 0) -
			2 * IM_GET_VOX(cur, x, y, z, 0) +
			IM_GET_VOX(prev, x, y, z, 0));

		// Solve the system
		switch(solve_Mat_rm(&Hs, &B, -1.0, &X)) {
		    case SUCCESS:
			break;
		    case SINGULAR:
			// The system is singular: give up
			goto refine_quit;
		    default:
			puts("refine_keypoint: error solving system! \n");
			return FAILURE;	
		}
	
#else
		// Parabolic interpolation
		MAT_RM_GET(&X, 0, 0, double) = -0.5 * ( 
		    IM_GET_VOX(cur, x + 1, y, z, 0) -
		    IM_GET_VOX(cur, x - 1, y, z, 0)) / (
		    IM_GET_VOX(cur, x + 1, y, z, 0) -
		    IM_GET_VOX(cur, x - 1, y, z, 0) +
		    2 * IM_GET_VOX(cur, x, y, z, 0));
		MAT_RM_GET(&X, 1, 0, double) = -0.5 * ( 
		    IM_GET_VOX(cur, x, y + 1, z, 0) -
		    IM_GET_VOX(cur, x, y - 1, z, 0)) / (
		    IM_GET_VOX(cur, x, y + 1, z, 0) -
		    IM_GET_VOX(cur, x, y - 1, z, 0) +
		    2 * IM_GET_VOX(cur, x, y, z, 0));
		MAT_RM_GET(&X, 2, 0, double) = -0.5 * ( 
		    IM_GET_VOX(cur, x, y, z + 1, 0) -
		    IM_GET_VOX(cur, x, y, z - 1, 0)) / (
		    IM_GET_VOX(cur, x, y, z + 1, 0) -
		    IM_GET_VOX(cur, x, y, z - 1, 0) +
		    2 * IM_GET_VOX(cur, x, y, z, 0));
		MAT_RM_GET(&X, 3, 0, double) = -0.5 * ( 
		    IM_GET_VOX(next, x, y, z, 0) -
		    IM_GET_VOX(prev, x, y, z, 0)) / (
		    IM_GET_VOX(next, x, y, z, 0) -
		    IM_GET_VOX(prev, x, y, z, 0) +
		    2 * IM_GET_VOX(cur, x, y, z, 0));
#endif
			// Update the coordinates
			xd = MAX(MIN(xd + MAT_RM_GET(&X, 0, 0, double), 
				     xmax), xmin);
			yd = MAX(MIN(yd + MAT_RM_GET(&X, 1, 0, double),	
				     ymax), ymin);
			zd = MAX(MIN(zd + MAT_RM_GET(&X, 2, 0, double), 
				     zmax), zmin);
			sd = MAX(MIN(sd + MAT_RM_GET(&X, 3, 0, double), 
				     smax), smin);

			// Compute the new pixel indices	
			xnew = (int) floor(xd);
			ynew = (int) floor(yd);
			znew = (int) floor(zd);

			// We are done if the pixel has not moved
			if (x == xnew && y == ynew && z == znew)
				break;

			// Update the pixel coordinates
			x = xnew;
			y = ynew;
			z = znew;
		}
#ifndef PARABOLA
refine_quit:
#else
#undef PARABOLA
#endif
		// Save the keypoint
		key->xi = x;
		key->yi = y;
		key->zi = z;
		key->xd = xd;
		key->yd = yd;
		key->zd = zd;
		key->sd = sd;
		key->sd_rel = sd * pow(2.0, -o);
	}

	return SUCCESS;
}

/* Bin a Cartesian gradient into Spherical gradient bins */
static int Cvec_to_sbins(const Cvec * const vd, Svec * const bins) {

	// Convert to spherical coordinates
	CVEC_TO_SVEC(vd, bins);
	//FIXME: Is this needed? CVEC_TO_SVEC cannot divide by zero
	if (bins->mag < FLT_EPSILON * 1E2)
		return FAILURE;

	// Compute bins
	bins->az *= (float) NBINS_AZ / AZ_MAX_F; 
	bins->po *= (float) NBINS_PO / PO_MAX_F;

	assert(bins->az < NBINS_AZ);
	assert(bins->po < NBINS_PO);

	return SUCCESS;
}

/* Refine an orientation histogram with optional operations,
 * such as solid angle weighting. */
static void refine_Hist(Hist *hist) {

#ifndef ICOS_HIST

#ifdef ORI_SOLID_ANGLE_WEIGHT
	{	
	float po;
	int a, p;
	// TODO: could accelerate this with a lookup table		

	// Weight by the solid angle of the bins, ignoring constants
	HIST_LOOP_START(a, p)
		po = p * (float) PO_MAX_F / NBINS_PO;
		HIST_GET(hist, a, p) /= cosf(po) - cosf(po + 
								(float) PO_MAX_F / NBINS_PO);
	HIST_LOOP_END
	}
#endif

#endif

}

/* Use orientation histograms to select all of the dominant
 * orientations for a keypoint. */
IGNORE_UNUSED
static int assign_hist_ori(const Image *const im, const Cvec *const vcenter,
                          const double sigma, Mat_rm *const R) {

	Svec buf[HIST_NUMEL];
	Hist hist, tmp_hist;
	Svec bins, temp;
	Cvec vd, vdisp, vx, vy, vz;
	Svec *ori;
	float sq_dist, weight, max, cur, af, pf, norm, proj;
	int i, x, y, z, a, p, n_ori;

	const double win_radius = sigma * ORI_RAD_FCTR;
	const float smooth_kernel[] = {1.0f / 16.0f, 4.0f / 16.0f, 
				       6.0f / 16.0f};

	// Initialize all magnitudes to -1
	for (i = 0; i < HIST_NUMEL; i++) {
		buf[i].mag = -1.0f;
	}

	// Initialize all histogram bins to 0
	HIST_LOOP_START(a, p)
		HIST_GET(&tmp_hist, a, p) = 0.0;
	HIST_LOOP_END

	IM_LOOP_SPHERE_START(im, x, y, z, vcenter, win_radius, &vdisp, sq_dist)

		// Compute Gaussian weighting, ignoring constant factor
		weight = expf(-0.5 * sq_dist / (sigma * sigma));    	
    
		// Get the gradient and convert to bins
		IM_GET_GRAD(im, x, y, z, 0, &vd);
		if (Cvec_to_sbins(&vd, &bins))
			continue;

		// Accumulate into hist
		a = (int) bins.az;
		p = (int) bins.po;
		HIST_GET(&tmp_hist, a, p) += bins.mag * weight; 	

		// TODO: Might benefit from trilinear interp. support	
			
	IM_LOOP_END

	// Histogram refinement steps
	refine_Hist(&tmp_hist);

	// Smooth the histogram in the azimuth angle (a la OpenCV)
	HIST_LOOP_START(a, p)
		HIST_GET(&hist, a, p) = (HIST_GET_AZ(&tmp_hist, a - 2, p) +
		   		HIST_GET_AZ(&tmp_hist, a + 2, p)) *	smooth_kernel[0] + 
	    	(HIST_GET_AZ(&tmp_hist, a - 1, p) +
		    	HIST_GET_AZ(&tmp_hist, a + 1, p)) * smooth_kernel[1] +
	    	HIST_GET(&tmp_hist, a, p) * smooth_kernel[2];
	HIST_LOOP_END

	// Smooth the histogram in the polar angle
	/*FIXME: This is not correct because HIST_GET_PO is not circular!!!
	  	 Smoothing in this way will bias results towards the endpoints 	*/
	HIST_LOOP_START(a, p)
		HIST_GET(&hist, a, p) = (HIST_GET_PO(&tmp_hist, a, p - 2) +
				HIST_GET_PO(&tmp_hist, a, p + 2)) * smooth_kernel[0] + 
			(HIST_GET_PO(&tmp_hist, a, p - 1) + 
				HIST_GET_PO(&tmp_hist, a, p + 1)) * smooth_kernel[1] +
			HIST_GET(&tmp_hist, a, p) * smooth_kernel[2];
	HIST_LOOP_END

	// Find the global maximum
	max = 0.0f;
	HIST_LOOP_START(a, p)
		max = MAX(HIST_GET(&hist, a, p), max);
	HIST_LOOP_END

	// Choose dominant orientations
	n_ori = 0;
	HIST_LOOP_START(a, p)

		// All orientations must be local extrema and within a
		// percentage of the global max  
		cur = HIST_GET(&hist, a, p);
		if (cur > ORI_LB * max && 
			cur > HIST_GET_AZ(&hist, a - 1, p    ) &&
			cur > HIST_GET_AZ(&hist, a - 1, p - 1) &&
			cur > HIST_GET_AZ(&hist, a - 1, p + 1) &&
			cur > HIST_GET_AZ(&hist, a + 1, p    ) &&
			cur > HIST_GET_AZ(&hist, a + 1, p - 1) &&
			cur > HIST_GET_AZ(&hist, a + 1, p + 1) &&
			cur > HIST_GET_PO(&hist, a,     p - 1) &&
			cur > HIST_GET_PO(&hist, a,     p + 1)) {

			af = (float) a;
			pf = (float) p;

#ifdef PARABOLOID_INTERP
			{
			float prev, next;

			// Interpolate the peak of an ellipic paraboloid.
			// See the appendix for more info. 
			next = HIST_GET_AZ(&hist, a + 1, p);
			prev = HIST_GET_AZ(&hist, a - 1, p);
			af += -0.5f * (next - prev) / (next - prev + 2.0f * cur);	
			af = fmodf(af + (float) NBINS_AZ, (float) NBINS_AZ);

			next = HIST_GET_PO(&hist, a, p + 1);
			prev = HIST_GET_PO(&hist, a, p - 1);
			pf += -0.5f * (next - prev) / (next - prev + 2.0f * cur);	
			pf = fmodf(pf + (float) NBINS_PO, (float) NBINS_PO);

			}
#endif

			// Save the dominant orientation
			ori = buf + n_ori++;
			ori->mag = cur;
			ori->az = af * AZ_MAX_F / (float) NBINS_AZ;
			ori->po = pf * PO_MAX_F / (float) NBINS_PO; 
			// TODO: Use quadratic interpolation for peak in po, az

			assert(ori->az < AZ_MAX_F);
			assert(ori->po < PO_MAX_F);
		}
	HIST_LOOP_END

	// Reject this keypoint if we have insufficient orientations 
	if (n_ori < 2) {
		return REJECT;
	}
	
	// Sort the orientations by decreasing magnitude
	qsort(buf, n_ori, sizeof(Svec), (__compar_fn_t) Svec_compar);

	// Test if the three strongest orientations are too similar
	if (buf[1].mag > ORI_UB * buf[0].mag || 
	    buf[2].mag > ORI_UB * buf[1].mag)
		return REJECT;
	
	// Convert orientations, and assign x' to the first one 
	temp.mag = 1.0f;
	temp.az = buf[0].az;
	temp.po = buf[0].po;
	SVEC_TO_CVEC(&temp, &vx);	
	temp.az = buf[1].az;
	temp.po = buf[1].po;
	SVEC_TO_CVEC(&temp, &vy);		

	// Perform Gram-Schmidt Orthogonalization to make y'
	proj = CVEC_DOT(&vx, &vy);
	vy.x -= proj * vx.x;
	vy.y -= proj * vx.y;
	vy.z -= proj * vx.z;
	norm = CVEC_L2_NORM(&vy);
	vy.x /= norm;
	vy.y /= norm;
	vy.z /= norm;

	// Compute z' = (x' X y')
	CVEC_CROSS(&vx, &vy, &vz);

	// Populate the rotation matrix
	MAT_RM_GET(R, 0, 0, float) = vx.x;
	MAT_RM_GET(R, 1, 0, float) = vx.y;
	MAT_RM_GET(R, 2, 0, float) = vx.z;

	MAT_RM_GET(R, 0, 1, float) = vy.x;
	MAT_RM_GET(R, 1, 1, float) = vy.y;
	MAT_RM_GET(R, 2, 1, float) = vy.z;

	MAT_RM_GET(R, 0, 2, float) = vz.x;
	MAT_RM_GET(R, 1, 2, float) = vz.y;
	MAT_RM_GET(R, 2, 2, float) = vz.z;

	return SUCCESS;
}

/* As above, but using the eigenvector method */
IGNORE_UNUSED
static int assign_eig_ori(const Image *const im, const Cvec *const vcenter,
                          const double sigma, Mat_rm *const R) {

    Cvec v[2];
    Mat_rm A, L, Q;
    Cvec vd, vd_win, vdisp, vr;
    double d, cos_ang;
    float weight, sq_dist, sgn;
    int i, x, y, z, m;
  
    const double win_radius = sigma * ORI_RAD_FCTR; 

    // Initialize the intermediates
    if (init_Mat_rm(&A, 3, 3, DOUBLE, TRUE) ||
	init_Mat_rm(&L, 0, 0, DOUBLE, TRUE) ||
	init_Mat_rm(&Q, 0, 0, DOUBLE, TRUE))
	goto eig_ori_fail;

    // Form the structure tensor and window gradient
    vd_win.x = 0.0f;
    vd_win.y = 0.0f;
    vd_win.z = 0.0f;
    IM_LOOP_SPHERE_START(im, x, y, z, vcenter, win_radius, &vdisp, sq_dist)
	// Compute Gaussian weighting, ignoring constant factor
	weight = expf(-0.5 * sq_dist / (sigma * sigma));		

	// Get the gradient	
	IM_GET_GRAD(im, x, y, z, 0, &vd);

	// Update the structure tensor
	MAT_RM_GET(&A, 0, 0, double) += (double) vd.x * vd.x * weight;
	MAT_RM_GET(&A, 0, 1, double) += (double) vd.x * vd.y * weight;
	MAT_RM_GET(&A, 0, 2, double) += (double) vd.x * vd.z * weight;
	MAT_RM_GET(&A, 1, 1, double) += (double) vd.y * vd.y * weight;
	MAT_RM_GET(&A, 1, 2, double) += (double) vd.y * vd.z * weight;
	MAT_RM_GET(&A, 2, 2, double) += (double) vd.z * vd.z * weight;

	// Update the window gradient
	CVEC_OP(&vd_win, &vd, +, &vd_win);
    IM_LOOP_END

    // Fill in the remaining elements
    MAT_RM_GET(&A, 1, 0, double) = MAT_RM_GET(&A, 0, 1, double);
    MAT_RM_GET(&A, 2, 0, double) = MAT_RM_GET(&A, 0, 2, double);
    MAT_RM_GET(&A, 2, 1, double) = MAT_RM_GET(&A, 1, 2, double);

    // Reject keypoints with weak gradient 
    if (CVEC_L2_NORM_SQ(&vd_win) < (float) ORI_GRAD_THRESH) {
	goto eig_ori_reject;
    } 

    // Get the eigendecomposition
    if (eigen_Mat_rm(&A, &Q, &L))
	goto eig_ori_fail;

    // Ensure we have distinct eigenvalues
    m = L.num_rows;
    if (m != 3)
	goto eig_ori_reject;

    // Test the corner response of the keypoint
    if (fabs(MAT_RM_GET(&L, 0, 0, double)) < EIG_MIN)
	goto eig_ori_reject;

    // Test the eigenvectors for stability
    for (i = 0; i < m - 1; i++) {
	if (fabs(MAT_RM_GET(&L, i, 0, double) /
		 MAT_RM_GET(&L, i + 1, 0, double)) > EIG_MAX_RATIO)
	    goto eig_ori_reject;
    }

    // Assign signs to the first n - 1 vectors
    for (i = 0; i < m - 1; i++) {

	const int eig_idx = m - i - 1;

	// Get an eigenvector, in descending order
	vr.x = (float) MAT_RM_GET(&Q, 0, eig_idx, double);
	vr.y = (float) MAT_RM_GET(&Q, 1, eig_idx, double);
	vr.z = (float) MAT_RM_GET(&Q, 2, eig_idx, double);

	// Get the directional derivative
	d = CVEC_DOT(&vd, &vr);

        // Get the cosine of the angle between the eigenvector and the gradient
        cos_ang = d / (CVEC_L2_NORM(&vr) * CVEC_L2_NORM(&vd));

        if (fabs(cos_ang) < EIG_COS_ANGLE_MIN) 
                goto eig_ori_reject;

	// Get the sign of the derivative
	if (d > 0.0)
	    sgn = 1.0f;
	else
	    sgn = -1.0f;

	// Enforce positive directional derivative
	CVEC_SCALE(&vr, sgn);

	// Add the vector to the rotation matrix
	MAT_RM_GET(R, 0, i, float) = vr.x;
	MAT_RM_GET(R, 1, i, float) = vr.y;
	MAT_RM_GET(R, 2, i, float) = vr.z;

	// Save this vector for later use
	v[i] = vr;
    }

    // Take the cross product of the first two vectors
    CVEC_CROSS(v, v + 1, &vr);

    // Add the last vector
    MAT_RM_GET(R, 0, 2, float) = (float) vr.x;
    MAT_RM_GET(R, 1, 2, float) = (float) vr.y;
    MAT_RM_GET(R, 2, 2, float) = (float) vr.z;

    cleanup_Mat_rm(&A);
    cleanup_Mat_rm(&Q);
    cleanup_Mat_rm(&L);
    return SUCCESS; 

eig_ori_reject:
    cleanup_Mat_rm(&A);
    cleanup_Mat_rm(&Q);
    cleanup_Mat_rm(&L);
    return REJECT;

eig_ori_fail:
    cleanup_Mat_rm(&A);
    cleanup_Mat_rm(&Q);
    cleanup_Mat_rm(&L);
    return FAILURE;
}

/* Assign rotation matrices to the keypoints. 
 * 
 * Note that this stage will modify kp, likely
 * rejecting some keypoints as orientationally
 * unstable. */
static int assign_orientations(SIFT3D *sift3d, 
			       Keypoint_store *kp) {

	int (*ori_fun)(const Image *const, const Cvec *const, const double, 
                       Mat_rm *const);
	Keypoint *kp_pos;
	size_t num;
	int i; 

    // Choose the orientation function to use
#ifdef EIG_ORI
    ori_fun = assign_eig_ori;
#else
    ori_fun = assign_hist_ori;
#endif

	// Iterate over the keypoints 
	kp_pos = kp->buf;
	for (i = 0; i < kp->slab.num; i++) {

		Keypoint *const key = kp->buf + i;
		const Image *const level = 
                        PYR_IM_GET(&sift3d->gpyr, key->o, key->s);
                Mat_rm *const R = &key->R;
                const Cvec vcenter = {key->xd, key->yd, key->zd};
                const double sigma = ORI_SIG_FCTR * key->sd_rel;

		// Initialize the orientation matrix
		if (init_Mat_rm_p(R, key->r_data, 3, 3, FLOAT, FALSE))
			return FAILURE;

		// Compute dominant orientations
		switch (ori_fun(level, &vcenter, sigma, R)) {
			case SUCCESS:
				// Continue processing this keypoint
				break;
			case REJECT:
				// Skip this keypoint
				continue;
			default:
				// Any other return value is an error
				return FAILURE;
		}
		
		// Rebuild the Keypoint buffer in place
		*kp_pos++ = *key; 
	}

	// Release unneeded keypoint memory
	num = kp_pos - kp->buf;
	RESIZE_KP_STORE(kp, num, sizeof(Keypoint));

	return SUCCESS;
}

/* Detect keypoint locations and orientations. You must initialize
 * the SIFT3D struct, image, and keypoint store with the appropriate
 * functions prior to calling this function. */
int SIFT3D_detect_keypoints(SIFT3D *const sift3d, Image *const im, 
			    Keypoint_store *const kp) {

	Image im_old; 

        // Verify inputs
        if (im->nc != 1) {
                fprintf(stderr, "SIFT3D_detect_keypoints: invalid number "
                        "of image channels: %d -- only single-channel images "
                        "are supported \n", im->nc);
                return FAILURE;
        }

	// Copy the last image and load the new one
	if (sift3d->im == NULL)
		// The SIFT3D has not been used: initialize im_old
		// to dummy values to mark for resizing.
		im_old.nx = im_old.ny = im_old.nz = -1;
	else
		im_old = *sift3d->im;
	sift3d->im = im;

	// Resize sift3d, if necessary
	if ((im->nx != im_old.nx || im->ny != im_old.ny ||
		im->nz != im_old.nz) && 
		resize_SIFT3D(sift3d))
		return FAILURE;

	// Build the GSS pyramid
	if (build_gpyr(sift3d))
		return FAILURE;

	// Build the DoG pyramid
	if (build_dog(sift3d))
		return FAILURE;

	// Detect extrema
	if (detect_extrema(sift3d, kp))
		return FAILURE;

	// Refine keypoints	
	if (refine_keypoints(sift3d, kp))
		return FAILURE;

	// Assign orientations
	if (assign_orientations(sift3d, kp))
		return FAILURE;

	return SUCCESS;
}

/* Get the bin and barycentric coordinates of a vector in the icosahedral 
 * histogram. */
IGNORE_UNUSED
static int icos_hist_bin(const SIFT3D * const sift3d,
			   const Cvec * const x, Cvec * const bary,
			   int * const bin) { 

	float k;
	int i;

	const Mesh * const mesh = &sift3d->mesh;

	// Check for very small vectors
	if (CVEC_L2_NORM_SQ(x) < BARY_EPS)
		return FAILURE;

	// Iterate through the faces
	for (i = 0; i < ICOS_NFACES; i++) {

		const Tri * const tri = mesh->tri + i;

		// Convert to barycentric coordinates
		if (cart2bary(x, tri, bary, &k))
			continue;

		// Test for intersection
		if (bary->x < -BARY_EPS || bary->y < -BARY_EPS || 
		    bary->z < -BARY_EPS || k < 0)
			continue;

		// Save the bin
		*bin = i;

		// No other triangles will be intersected
		return SUCCESS;
	}	

	// Unreachable code
	assert(FALSE);
	return FAILURE;
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
	mag = CVEC_L2_NORM(grad);

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
				weight = ((dx == 0) ? (1.0f - dvbins.x) : 
					 	dvbins.x) *
					 ((dy == 0) ? (1.0f - dvbins.y) : 
						dvbins.y) *
					 ((dz == 0) ? (1.0f - dvbins.z) : 
						dvbins.z);

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
					}
				}	
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
static void extract_descrip(SIFT3D *sift3d, Image *im,
	   Keypoint *key, SIFT3D_Descriptor *desc) {

	Cvec vcenter, vim, vkp, vbins, grad, grad_rot;
	Hist *hist;
	float weight, sq_dist;
	int i, x, y, z, a, p;

	// Compute basic parameters 
        const float sigma = key->sd_rel;
	const float win_radius = DESC_RAD_FCTR * sigma;
	const float desc_width = win_radius / SQRT2_F;
	const float desc_hw = desc_width / 2.0f;
	const float desc_bin_fctr = (float) NHIST_PER_DIM / desc_width;
	const double coord_factor = pow(2.0, key->o);

	// Zero the descriptor
	for (i = 0; i < DESC_NUM_TOTAL_HIST; i++) {
		hist = desc->hists + i;
                hist_zero(hist);
	}

	// Iterate over a sphere window in image space
	vcenter.x = key->xd;
	vcenter.y = key->yd;
	vcenter.z = key->zd;
	IM_LOOP_SPHERE_START(im, x, y, z, &vcenter, win_radius, &vim, sq_dist)

		// Rotate to keypoint space
		MUL_MAT_RM_CVEC(&key->R, &vim, &vkp);		

		// Compute spatial bins
		vbins.x = (vkp.x + desc_hw) * desc_bin_fctr;
		vbins.y = (vkp.y + desc_hw) * desc_bin_fctr;
		vbins.z = (vkp.z + desc_hw) * desc_bin_fctr;

		// Reject points outside the rectangular descriptor 
		if (vbins.x < 0 || vbins.y < 0 || vbins.z < 0 ||
			vbins.x >= (float) NHIST_PER_DIM ||
			vbins.y >= (float) NHIST_PER_DIM ||
			vbins.z >= (float) NHIST_PER_DIM)
			continue;

		// Take the gradient
		IM_GET_GRAD(im, x, y, z, 0, &grad);

		// Apply a Gaussian window
		weight = expf(-0.5f * sq_dist / (sigma * sigma));
		CVEC_SCALE(&grad, weight);

                // Rotate the gradient to keypoint space
		MUL_MAT_RM_CVEC(&key->R, &grad, &grad_rot);

		// Finally, accumulate bins by 5x linear interpolation
		SIFT3D_desc_acc_interp(sift3d, &vbins, &grad_rot, desc);
	IM_LOOP_END

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
			HIST_GET(hist, a, p) = MIN(HIST_GET(hist, a, p), 
						   TRUNC_THRESH);
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
}

/* Extract SIFT3D descriptors from a list of keypoints and an 
 * image.
 *
 * parameters:
 * 	sift3d - (initialized) struct defining parameters
 *  im - If use_gpyr == 0, this is a pointer to an Image
 *   	 and features will be extracted from the "raw" data.
 *  	 Else, this is a pointer to a Pyramid and features
 * 		 will be extracted at the nearest pyramid level to
 * 		 the keypoint. 
 *  kp - keypoint list populated by a feature detector 
 *  desc - (initialized) struct to hold the descriptors
 *  use_gpyr - see im for details */
int SIFT3D_extract_descriptors(SIFT3D *sift3d, void *im,
	Keypoint_store *kp,	SIFT3D_Descriptor_store *desc, int use_gpyr) {

	Image *level;
	SIFT3D_Descriptor *descrip;
	Keypoint *key;
	int i;

	// Parse inputs
	if (!use_gpyr) {
		puts("SIFT3D_extract_descriptors: This feature has not yet "
			 "been implemented! Please call this function with a "
			 "Pyramid rather than an Image. \n");
	}

	// Resize the descriptor store
	desc->num = kp->slab.num;
	if ((desc->buf = (SIFT3D_Descriptor *) realloc(desc->buf, desc->num * 
				sizeof(SIFT3D_Descriptor))) == NULL)
		return FAILURE;

	// Initialize the image info
	if (use_gpyr) {
		Pyramid *gpyr;
		gpyr = (Pyramid *) im;
		level = PYR_IM_GET(gpyr, gpyr->first_octave, gpyr->first_level);
	} else
		level = im;
	desc->nx = level->nx;	
	desc->ny = level->ny;	
	desc->nz = level->nz;	

        // Extract the descriptors
	for (i = 0; i < desc->num; i++) {
		key = kp->buf + i;
		descrip = desc->buf + i;
		level = PYR_IM_GET((Pyramid *) im, key->o, key->s);
		extract_descrip(sift3d, level, key, descrip);
	}	

	return SUCCESS;
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
static void postproc_Hist(Hist *const hist) {

        int a, p;

        const float hist_trunc = TRUNC_THRESH * DESC_NUMEL / HIST_NUMEL;

	// Histogram refinement steps
	refine_Hist(hist);

	// Normalize the descriptor
	normalize_hist(hist);

	// Truncate
	HIST_LOOP_START(a, p)
		HIST_GET(hist, a, p) = MIN(HIST_GET(hist, a, p), 
					   hist_trunc);
	HIST_LOOP_END

	// Normalize again
	normalize_hist(hist);
}

/* Helper routine to extract a single SIFT3D histogram. */
static void extract_dense_descrip(SIFT3D *const sift3d, 
           const Image *const im, const Cvec *const vcenter, 
           const double sigma, const Mat_rm *const R, Hist *const hist) {

	Cvec grad, grad_rot, bary, vim;
	float sq_dist, mag, weight;
	int a, p, x, y, z, bin;

        const Mesh *const mesh = &sift3d->mesh;
	const float win_radius = DESC_RAD_FCTR * sigma;
	const float desc_width = win_radius / SQRT2_F;
	const float desc_hw = desc_width / 2.0f;

	// Zero the descriptor
        hist_zero(hist);

	// Iterate over a sphere window in image space
	IM_LOOP_SPHERE_START(im, x, y, z, vcenter, win_radius, &vim, sq_dist)

		// Take the gradient and rotate
		IM_GET_GRAD(im, x, y, z, 0, &grad);
		MUL_MAT_RM_CVEC(R, &grad, &grad_rot);

                // Get the index of the intersecting face
                if (icos_hist_bin(sift3d, &grad_rot, &bary, &bin))
                        continue;

                // Get the magnitude of the vector
                mag = CVEC_L2_NORM(&grad);

		// Get the Gaussian window weight
		weight = expf(-0.5f * sq_dist / (sigma * sigma));

                // Interpolate over three vertices
                MESH_HIST_GET(mesh, hist, bin, 0) += mag * weight * bary.x;
                MESH_HIST_GET(mesh, hist, bin, 1) += mag * weight * bary.y;
                MESH_HIST_GET(mesh, hist, bin, 2) += mag * weight * bary.z;

	IM_LOOP_END

        // Histogram postprocessing
        postproc_Hist(hist);
}

/* Get a descriptor with a single histogram at each voxel of an image.
 * The result is an image with HIST_NUMEL channels, where each channel is a
 * bin of the histogram. 
 * Parameters:
 * -sift3d The descriptor extractor.
 * -in The input image.
 * -out The output image.
 */
int SIFT3D_extract_dense_descriptors(SIFT3D *const sift3d, 
        const Image *const in, Image *const desc) {

        // Verify inputs
        if (in->nc != 1) {
                fprintf(stderr, "SIFT3D_extract_dense_descriptors: invalid "
                        "number of channels: %d. This function only supports "
                        "single-channel images. \n", in->nc);
                return FAILURE;
        }

        // Resize the output image
        memcpy(desc->dims, in->dims, IM_NDIMS * sizeof(int));
        desc->nc = HIST_NUMEL;
        im_default_stride(desc);
        if (im_resize(desc))
                return FAILURE;

        // Extract the descriptors
        return sift3d->dense_rotate ? 
                extract_dense_descriptors_rotate(sift3d, in, desc) :
                extract_dense_descriptors(sift3d, in, desc);
}

/* Helper function for extract_dense_descriptors, without rotation invariance.
 * This function is much faster than its rotation-invariant counterpart 
 * because histogram bins are pre-computed. */
static int extract_dense_descriptors(SIFT3D *const sift3d,
        const Image *const in, Image *const desc) {

        Image temp; 
        Gauss_filter gauss;
	Cvec grad, grad_rot, bary;
	float sq_dist, mag, weight;
        int i, x, y, z, bin, vert;

        const int x_start = 1;
        const int y_start = 1;
        const int z_start = 1;
        const int x_end = in->nx - 2;
        const int y_end = in->ny - 2;
        const int z_end = in->nz - 2;

        Mesh * const mesh = &sift3d->mesh;
        const double sigma_win = sift3d->dense_sigma * DESC_SIG_FCTR;

        // Initialize the intermediate image
        init_im(&temp);
        if (im_copy_dims(desc, &temp))
                return FAILURE;

        // Initialize the filter
        if (init_Gauss_filter(&gauss, sigma_win, 3)) {
                im_free(&temp);
                return FAILURE;
        }

        // Initialize the descriptors for each voxel
        im_zero(&temp);
        IM_LOOP_LIMITED_START(in, x, y, z, x_start, x_end, y_start, y_end, 
                              z_start, z_end)

                // Take the gradient
		IM_GET_GRAD(in, x, y, z, 0, &grad);

                // Get the index of the intersecting face
                if (icos_hist_bin(sift3d, &grad_rot, &bary, &bin))
                        continue;

                // Initialize each vertex
                IM_GET_VOX(&temp, x, y, z, MESH_GET_IDX(mesh, bin, 0)) = bary.x;
                IM_GET_VOX(&temp, x, y, z, MESH_GET_IDX(mesh, bin, 1)) = bary.y;
                IM_GET_VOX(&temp, x, y, z, MESH_GET_IDX(mesh, bin, 2)) = bary.z;

        IM_LOOP_END

        // Filter the descriptors
	if (apply_Sep_FIR_filter(&temp, desc, &gauss.f))
                goto dense_extract_quit;

        // Post-process the descriptors
        IM_LOOP_START(desc, x, y, z)

                Hist hist;

                // Copy to a Hist
                vox2hist(desc, x, y, z, &hist);

                // Post-process
                postproc_Hist(&hist);

                // Copy back to the image
                hist2vox(&hist, desc, x, y, z);

        IM_LOOP_END

        // Clean up
        im_free(&temp);
        cleanup_Gauss_filter(&gauss);

        return SUCCESS;

dense_extract_quit:
        im_free(&temp);
        cleanup_Gauss_filter(&gauss);
        return FAILURE;
}

/* Copy a voxel to a Hist. Does no bounds checking. */
static int vox2hist(const Image *const im, const int x, const int y,
        const int z, Hist *const hist) {

        int c;

        for (c = 0; c < HIST_NUMEL; c++) {
                hist->bins[c] = IM_GET_VOX(im, x, y, z, c);
        }
}

/* Copy a Hist to a voxel. Does no bounds checking. */
static int hist2vox(Hist *const hist, const Image *const im, const int x, 
        const int y, const int z) {

        int c;
        
        for (c = 0; c < HIST_NUMEL; c++) {
                IM_GET_VOX(im, x, y, z, c) = hist->bins[c];
        }
}

/* As in extract_dense_descrip, but with rotation invariance */
static int extract_dense_descriptors_rotate(SIFT3D *const sift3d,
        const Image *const in, Image *const desc) {

        Hist hist;
        Mat_rm R, Id;
        Mat_rm *ori;
        int i, x, y, z, c;

        // Initialize the identity matrix
        if (init_Mat_rm(&Id, 3, 3, FLOAT, TRUE)) {
                return FAILURE;
        }
        for (i = 0; i < 3; i++) {       
                MAT_RM_GET(&Id, i, i, float) = 1.0f;
        }

        // Initialize the rotation matrix
        if (init_Mat_rm(&R, 3, 3, FLOAT, TRUE)) {
                cleanup_Mat_rm(&Id);
                return FAILURE;
        }

        // Iterate over each voxel
        IM_LOOP_START(in, x, y, z)

                const Cvec vcenter = {(float) x + 0.5f, 
                                      (float) y + 0.5f, 
                                      (float) z + 0.5f};

                const double sigma = sift3d->dense_sigma;

                // Attempt to assign an orientation
                switch(assign_eig_ori(in, &vcenter, sigma, &R)) {
                        case SUCCESS:
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
                extract_dense_descrip(sift3d, in, &vcenter, sigma, ori,
                                      &hist);

                // Copy the descriptor to the image channels
                hist2vox(&hist, desc, x, y, z);

        IM_LOOP_END

        // Clean up
        cleanup_Mat_rm(&R);
        cleanup_Mat_rm(&Id);
        return SUCCESS;

dense_rotate_quit:
        // Clean up and return an error condition 
        cleanup_Mat_rm(&R);
        cleanup_Mat_rm(&Id);
        return FAILURE;
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
    mat->type = DOUBLE;
    if (resize_Mat_rm(mat))
	return FAILURE;

    // Build the matrix
    for (i = 0; i < num; i++) {
	MAT_RM_GET(mat, i, 0, double) = kp->buf[i].xd;
	MAT_RM_GET(mat, i, 1, double) = kp->buf[i].yd;
	MAT_RM_GET(mat, i, 2, double) = kp->buf[i].zd;
    }

    return SUCCESS;
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
		return FAILURE;
	}

	// Resize inputs
	mat->type = FLOAT;
	mat->num_rows = num_rows;
	mat->num_cols = num_cols;
	if (resize_Mat_rm(mat))
		return FAILURE;

	// Copy the data
	for (i = 0; i < num_rows; i++) {

		const SIFT3D_Descriptor *const desc = store->buf + i;

		// Copy the coordinates
		MAT_RM_GET(mat, i, 0, float) = (float) desc->xd;
		MAT_RM_GET(mat, i, 1, float) = (float) desc->yd;
		MAT_RM_GET(mat, i, 2, float) = (float) desc->zd;

		// Copy the feature vector
		for (j = 0; j < DESC_NUM_TOTAL_HIST; j++) {
			const Hist *const hist = desc->hists + j;
			HIST_LOOP_START(a, p)
				MAT_RM_GET(mat, i, j + IM_NDIMS, float) = 
					HIST_GET(hist, a, p);
			HIST_LOOP_END
		}
	}

	return SUCCESS;
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
		printf("Mat_rm_to_SIFT3D_Descriptor_store: invalid matrix "
		       "dimensions: [%d X %d] \n", num_rows, num_cols);
		return FAILURE;
	}

	// Resize the descriptor store
	store->num = num_rows;
	if ((store->buf = (SIFT3D_Descriptor *) realloc(store->buf, store->num * 
				sizeof(SIFT3D_Descriptor))) == NULL)
		return FAILURE;

	// Copy the data
	for (i = 0; i < num_rows; i++) {

		SIFT3D_Descriptor *const desc = store->buf + i;

		// Copy the coordinates
		desc->xd = MAT_RM_GET(mat, i, 0, float);
		desc->yd = MAT_RM_GET(mat, i, 1, float);
		desc->zd = MAT_RM_GET(mat, i, 2, float);

		// Copy the feature vector
		for (j = 0; j < DESC_NUM_TOTAL_HIST; j++) {
			Hist *const hist = desc->hists + j;
			HIST_LOOP_START(a, p)
				HIST_GET(hist, a, p) = MAT_RM_GET(mat, i, 
					j + IM_NDIMS, float);
			HIST_LOOP_END
		}
	}
	
	return SUCCESS;
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
    match1->type = match2->type = DOUBLE;
    if (resize_Mat_rm(match1) || resize_Mat_rm(match2))
	    return FAILURE;

    // Populate the matrices
    num_matches = 0;
    for (i = 0; i < num; i++) {

	    const SIFT3D_Descriptor *const desc1 = d1->buf + i;
	    const SIFT3D_Descriptor *const desc2 = d2->buf + matches[i];

	    if (matches[i] == -1)
		    continue;

	    // Save the match
	    MAT_RM_GET(match1, num_matches, 0, double) = desc1->xd; 
	    MAT_RM_GET(match1, num_matches, 1, double) = desc1->yd; 
	    MAT_RM_GET(match1, num_matches, 2, double) = desc1->zd; 
	    MAT_RM_GET(match2, num_matches, 0, double) = desc2->xd; 
	    MAT_RM_GET(match2, num_matches, 1, double) = desc2->yd; 
	    MAT_RM_GET(match2, num_matches, 2, double) = desc2->zd; 
	    num_matches++;
    }

    // Release extra memory
    match1->num_rows = match2->num_rows = num_matches;
    if (resize_Mat_rm(match1) || resize_Mat_rm(match2))
	    return FAILURE;
    
    return SUCCESS;
}

/* Like SIFT3D_nn_match, but also tests for forward-backward consistency.
 * That is, matching is performed from d1 to d2, and then d2 to d1, and
 * any matches that do not have the same result in each pass are rejected. 
 *
 * See SIFT3D_nn_match for descriptions of the parameters. */
int SIFT3D_nn_match_fb(const SIFT3D_Descriptor_store *const d1,
		       const SIFT3D_Descriptor_store *const d2,
		       const float nn_thresh, int **const matches) {
    int *matches2; 
    int i;

    // Run the matching passes
    matches2 = NULL;
    if (SIFT3D_nn_match(d1, d2, nn_thresh, matches) ||
	SIFT3D_nn_match(d2, d1, nn_thresh, &matches2))
	return FAILURE;

    // Enforce forward-backward consistency
    for (i = 0; i < d1->num; i++) {

	int *const match1 = *matches + i;

	if (*match1 >= 0 && matches2[*match1] != i)
	    *match1 = -1;
    }

    return SUCCESS;
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
 * You may consider using SIFT3D_matches_to_Mat_rm to convert the matches to
 * coordinate matrices. */
int SIFT3D_nn_match(const SIFT3D_Descriptor_store *const d1,
		    const SIFT3D_Descriptor_store *const d2,
		    const float nn_thresh, int **const matches) {

	const SIFT3D_Descriptor *desc_best;
	float *match_ssds;
	float ssd, ssd_best, ssd_nearest;
	int i, j, k, a, p, desc2_idx;

	const int num = d1->num;

#ifdef MATCH_MAX_DIST
		Cvec dims, dmatch;
		double dist_match;
				
		// Compute spatial distance rejection threshold
		dims.x = (float) d1->nx;	
		dims.y = (float) d1->ny;	
		dims.z = (float) d1->nz;	
		const double diag = CVEC_L2_NORM(&dims);	
		const double dist_thresh = diag * MATCH_MAX_DIST;
#endif

	// Initialize intermediate arrays 
	if ((*matches = (int *) realloc(*matches, num * sizeof(float))) == NULL ||
	    (match_ssds = (float *) malloc(num * sizeof(float))) == NULL) {
	    fprintf(stderr, "SIFT3D_nn_match: out of memory! \n");
	    return FAILURE;
	}

	for (i = 0; i < d1->num; i++) {
	    // Mark -1 to signal there is no match
	    (*matches)[i] = -1;
	    match_ssds[i] = -1.0f;
	}
	
	// Exhaustive search for matches
	for (i = 0; i < d1->num; i++) {

	    const SIFT3D_Descriptor *const desc1 = d1->buf + i;

	    // Linear search for the best and second-best SSD matches 
	    ssd_best = ssd_nearest = 1e30f;
	    desc_best = NULL;
	    for (j = 0; j < d2->num; j++) { 

		    const SIFT3D_Descriptor *const desc2 = d2->buf + j;

		    // Find the SSD of the two descriptors
		    ssd = 0.0f;
		    for (k = 0; k < DESC_NUM_TOTAL_HIST; k++) {

			    const Hist *const hist1 = &desc1->hists[k];
			    const Hist *const hist2 = &desc2->hists[k];

			    HIST_LOOP_START(a, p)
				    const float diff = HIST_GET(hist1, a, p) - 
					       HIST_GET(hist2, a, p); 
				    ssd += diff * diff; 		
			    HIST_LOOP_END
		    }

		    // Compare to best matches
		    if (ssd < ssd_best) {
			    desc_best = desc2; 
			    ssd_nearest = ssd_best;
			    ssd_best = ssd;
		    } else 
			    ssd_nearest = MIN(ssd_nearest, ssd);
	    }

	    // Reject match if nearest neighbor is too close 
	    if (ssd_best / ssd_nearest > nn_thresh * nn_thresh)
		    goto match_reject;

	    desc2_idx = desc_best - d2->buf;

#ifdef MATCH_MAX_DIST
	    // Compute the spatial distance of the match
	    dmatch.x = (float) desc_best->xd - desc1->xd; 
	    dmatch.y = (float) desc_best->yd - desc1->yd; 
	    dmatch.z = (float) desc_best->zd - desc1->zd; 
	    dist_match = (double) CVEC_L2_NORM(&dmatch);

	    // Reject matches of great distance
	    if (dist_match > dist_thresh)
		    goto match_reject;
#endif
			
	    // Save the match
	    (*matches)[i] = desc2_idx;
	    match_ssds[i] = ssd_best;

match_reject: ;
	}

	return SUCCESS;
}

int draw_matches(const Image *const src, const Image *const ref,
		 const Mat_rm *const match_src, const Mat_rm *const match_ref,
		 Image *const background, Image *const overlay) {

	Mat_rm match_ref_draw;
	int i;

	// Initialize intermediates		
	if (init_Mat_rm(&match_ref_draw, 0, 0, DOUBLE, FALSE) ||
	    copy_Mat_rm(match_ref, &match_ref_draw))
	return FAILURE;

	// Adjust the coordinates of the reference features 
	for (i = 0; i < match_ref->num_rows; i++) {
		MAT_RM_GET(&match_ref_draw, i, 0, double) =
			MAT_RM_GET(match_ref, i, 0, double) + (double) src->nx;
	}

	// Draw a concatenated image
	if (im_concat(src, ref, 0, background))
		return FAILURE;

	// Draw the matches
	if (draw_lines(match_src, &match_ref_draw, background->nx, 
		       background->ny, background->nz, overlay))
		return FAILURE;

	return SUCCESS;
}
