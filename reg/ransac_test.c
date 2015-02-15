/* ransac_test.c
 * ----------------------------------------------------------------
 * Rice MRI Team
 * ----------------------------------------------------------------
 * Unit test to get RANSAC affine registration working.
 * ----------------------------------------------------------------
 * Created: Daniel Reiter
 * Last updated: Daniel Reiter
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include "types.h"
#include "sift.h"
#include "macros.h"
#include "imutil.h"
#include "reg.h"

int main(int argc, char *argv[]) {

	Mat_rm match_src, match_ref, A;
	Affine aff_reg;
	Ransac ran;

	// Initialization
	if (init_Mat_rm(&match_src, 4, 3, DOUBLE, FALSE) ||
		init_Mat_rm(&match_ref, 4, 3, DOUBLE, FALSE) ||
		init_Mat_rm(&A, 4, 4, DOUBLE, TRUE))
		return FAILURE;
	init_Affine(&aff_reg, 3);
	init_Ransac(&ran, 0.3, 400, 20);

	// Set up affine transformation matrix to later compare with
	double ang_deg = 5.0;
	double ang_rad = ang_deg * UTIL_PI / 180.0;
	MAT_RM_GET(&A, 0, 0, double) = cos(ang_rad);
	MAT_RM_GET(&A, 0, 1, double) = -sin(ang_rad);
	MAT_RM_GET(&A, 1, 0, double) = sin(ang_rad);
	MAT_RM_GET(&A, 1, 1, double) = cos(ang_rad);
	MAT_RM_GET(&A, 2, 2, double) = 1.0;
	
	// initialize match_src and match_ref to simple precomputed transformation
	int num_pts = 4;
	int dim = 3;
	double match_ref_vals[] = {150, 50, 50,
				50, 150, 50,
				50, 50, 150,
				//25, 30, 249,
				//50, 15, 109,
				// TODO: need a 5th point from this transformation
				50, 150, 130};
	// 5 deg rotation
	double match_src_vals[] = {145.0714, 62.8831, 50.0000,
				36.7364, 153.7870, 50.0000,
				45.4519, 54.1675, 150.0000,
                                //1, 1, 1,
				//200, 105, 30
				36.7364, 153.7870, 130.0000};
	int i, j, index;
	for (i = 0; i < num_pts; i++) {
	    for (j = 0; j < dim; j++) {
		index = i*dim + j;
		MAT_RM_GET(&match_ref, i, j, double) = match_ref_vals[index];
		MAT_RM_GET(&match_src, i, j, double) = match_src_vals[index];
	    }
	}

	puts("Initial transformation: \n");
	print_Mat_rm(&A);

	// Find the transformation 
	printf("Computing transformation...\n");
	if (find_tform_ransac(&ran, &match_src, &match_ref, 3, AFFINE, 
						  (void *) &aff_reg))	
		err_exit("No good model found\n");

	puts("Inferred transformation: \n");
	print_Mat_rm(&aff_reg.A);

	return 0;
}
