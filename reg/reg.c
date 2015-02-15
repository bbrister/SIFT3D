/* reg.c
 * ----------------------------------------------------------------
 * Rice MRI Team
 * ----------------------------------------------------------------
 * This file contains all routines needed to do a RANSAC-based
 * affine image registration.
 *-----------------------------------------------------------------
 * Created: Jenny Zhang 1/24/2014
 * Last updated: Blaine Rister 11/18/2014
 */

#include "reg.h"
#include "imutil.h"

/* Helper routines and data structures. */
typedef struct _List {
  struct _List *next;
  struct _List *prev;
  int idx;
} List;

/* Unfinished public routines */
int init_Tps(Tps *tps, int dim, int terms);
int resize_Tps(Tps* tps, int num_pts, int dim);

/* Helper routines */
static int init_List(List **list, const int num);
static int List_get(List *list, const int idx, List **el);
static void List_remove(List **list, List *el);
static void cleanup_List(List *list);
static int rand_rows(Mat_rm *in1, Mat_rm *in2, Mat_rm *out1, Mat_rm *out2, 
              int num_rows);
static int make_spline_matrix(Mat_rm* src, Mat_rm* src_in, Mat_rm* sp_src, int K_terms, int* r, int dim);
static int make_affine_matrix(Mat_rm* pts_in, Mat_rm* mat_out, const int dim);
static Mat_rm* extract_ctrl_pts(void* tform, tform_type type);
static Mat_rm* extract_ctrl_pts_Tps(Tps *tps);
static int solve_system(void* tform, Mat_rm* src, Mat_rm* ref, const int dim, 
                 tform_type type);
static double tform_err_sq(void* tform, Mat_rm* src, Mat_rm* ref, int i, 
			   tform_type type);
static int ransac(Mat_rm* src, Mat_rm* ref, Ransac* ran, const int dim, void *tform, 
           tform_type type, int **cset, int *len);

/* Initialize an Tps struct. This initializes
 * all fields, and allocates memory for the inner
 * matrix, initializing it to zero. */
int init_Tps(Tps *tps, int dim, int terms) {
	// Verify inputs
	if (dim < 2)
		return FAILURE;
	
	// Initialize the matrix
	if (init_Mat_rm(&tps->params, dim, terms,
					DOUBLE, TRUE))
		return FAILURE;
		
	if (init_Mat_rm(&tps->kp_src, terms-dim-1,dim,
					DOUBLE, TRUE))
		return FAILURE;
	return SUCCESS;
	
	tps->dim=dim;
}

/* Initialize a RANSAC struct with the given parameters */
int init_Ransac(Ransac *ran, double min_inliers, double err_thresh, 
				 int num_iter) {
	ran->min_inliers = min_inliers; //ratio of inliers threshold for RANSAC concensus sets
 	ran->err_thresh = err_thresh; //error threshold for RANSAC inliers
	ran->num_iter = num_iter; //number of RANSAC iterations

	return SUCCESS;
}

/* Select a random subset of rows, length "num_rows".
 * This function resizes out. 
 * 
 * Returns an error if in->num_rows < num_rows.
 * Both input matrices must have type double and the same dimensions.
 *
 * All matrices must be initialized prior to calling 
 * this function.*/
static int rand_rows(Mat_rm *in1, Mat_rm *in2, Mat_rm *out1, Mat_rm *out2, 
              int num_rows) {

        List *row_indices, *el;
        int i, j, list_size, idx;

	const int num_rows_in = in1->num_rows;
        const int num_cols = in1->num_cols;
        const int num_remove = in1->num_rows - num_rows;

        // Verify inputs
        if (in2->num_rows != num_rows_in || 
            in2->num_cols != num_cols) {
          puts("rand_rows; inputs must have the same dimension \n");
          return FAILURE;
        }
        if (in1->num_rows > RAND_MAX) {
          puts("rand_rows: input matrix is too large \n");
          return FAILURE;
        }
        if (num_remove < 0) {
          puts("rand_rows: not enough rows in the matrix \n");
          return FAILURE;
        }
        if (in1->type != DOUBLE || in2->type != DOUBLE) {
          puts("rand_rows: inputs must have type int \n");
          return FAILURE;
        }

        // Resize the outputs
        out1->type = out2->type = in1->type;
        out1->num_rows = out2->num_rows = num_rows;
        out1->num_cols = out2->num_cols = num_cols;
        if (resize_Mat_rm(out1) || resize_Mat_rm(out2))
          return FAILURE;

        // Initialize a list with all of the row indices
        if (init_List(&row_indices, num_rows_in))
          return FAILURE;

        for (i = 0; i < num_rows_in; i++) {
          if (List_get(row_indices, i, &el))
            return FAILURE;
          el->idx = i;
        }

        // Remove random rows
        list_size = num_rows_in;
        for (i = 0; i < num_remove; i++) {
          // Draw a random number
          idx = rand() % list_size;

	  // Remove that element
          if (List_get(row_indices, idx, &el))
		return FAILURE;
          List_remove(&row_indices, el);
          list_size--;
        }

        // Build the output matrices
        el = row_indices;
	MAT_RM_LOOP_START(out1, i, j)
            MAT_RM_GET(out1, i, j, double) = MAT_RM_GET(in1, el->idx, j, double);
            MAT_RM_GET(out2, i, j, double) = MAT_RM_GET(in2, el->idx, j, double);
	    MAT_RM_LOOP_COL_END

          // Get the next row
          el = el->next;
	MAT_RM_LOOP_ROW_END

        // Clean up
        cleanup_List(row_indices);

	return SUCCESS;
}

/* Initialize a list of num elements. */
static int init_List(List **list, const int num) {

  List *prev, *cur;
  int i;

  prev = NULL;
  for (i = 0; i < num; i++) {
    if ((cur = (List *) malloc(sizeof(List))) == NULL)
      return FAILURE;

    if (i == 0)
      *list = cur;

    cur->next = NULL;
    cur->prev = prev;

    if (prev != NULL)
      prev->next = cur;

    prev = cur;
  }

  return SUCCESS;
}

/* Get the element idx in list. Returns a pointer to the element in el.
 * Returns FAILURE if NULL is reached before idx elements have been 
 * traversed. */
static int List_get(List *list, const int idx, List **el) {

  int i;

    // Verify inputs
    if (idx < 0)
	return FAILURE;

  // Traverse the list 
  for (i = 0; i < idx; i++) {
    list = list->next;

    if (list == NULL)
      return FAILURE;
  }

  *el = list;
  return SUCCESS;
}

/* Remove an element from its list, freeing its memory. */
static void List_remove(List **list, List *el) {
  if (list == NULL || el == NULL)
    return;

  if (el->next != NULL)
    el->next->prev = el->prev;

  if (el->prev != NULL)
    el->prev->next = el->next;

  if (*list == el)
	*list = el->next;

  free(el);
}

/* Free all of the memory for a list. */
static void cleanup_List(List *list) {

  while (1) {
    List_remove(&list, list->prev);

    if (list->next == NULL) {
      List_remove(&list, list);
      break;
    }

    list = list->next;
  }

}

//make the system matrix for spline
static int make_spline_matrix(Mat_rm* src, Mat_rm* src_in, Mat_rm* sp_src, int K_terms, int* r, int dim){
	int i,d;
	double x,y,z,x2,y2,z2,r_sq,U;
	src_in->type = DOUBLE;
	sp_src->type = DOUBLE;
 	if (init_Mat_rm(src_in, K_terms+dim+1, K_terms+dim+1, DOUBLE, TRUE)){
		return FAILURE;
	}
 	if (init_Mat_rm(sp_src, K_terms, dim, DOUBLE, TRUE)){
		return FAILURE;
	}	
 	for (i=0;i<K_terms;i++){
		//get the coordinate of current point
		switch (dim) {
			case 2:
				x = MAT_RM_GET(src, r[i], 0, double);
				y = MAT_RM_GET(src, r[i], 1, double);
				break;
			case 3:
				x = MAT_RM_GET(src, r[i], 0, double);
				y = MAT_RM_GET(src, r[i], 1, double);
				z = MAT_RM_GET(src, r[i], 2, double);
				break;
		}
		for (d=0;d<i;d++){
			//compute r
			switch (dim) {
				case 2:
					x2 = MAT_RM_GET(src, r[d], 0, double);
					y2 = MAT_RM_GET(src, r[d], 1, double);
					r_sq=(x-x2)*(x-x2)+(y-y2)*(y-y2);
					break;
				case 3:
					x2 = MAT_RM_GET(src, r[d], 0, double);
					y2 = MAT_RM_GET(src, r[d], 1, double);
					z2 = MAT_RM_GET(src, r[d], 2, double);
					r_sq=(x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2);
					break;
			}
			//compute U
			U=r_sq*log(r_sq);
			//construct K
			MAT_RM_GET(src_in, i, d, double)=U;
			MAT_RM_GET(src_in, d, i, double)=U;				
		}
		MAT_RM_GET(src_in, i, i, double)=0.0;
		//construct P and P'
		MAT_RM_GET(src_in, i, K_terms, double)=1.0;
		MAT_RM_GET(src_in, K_terms, i, double)=1.0;
		switch (dim) {
			case 2:
				MAT_RM_GET(src_in, i, K_terms+1, double)=x;
				MAT_RM_GET(src_in, i, K_terms+2, double)=y;
				MAT_RM_GET(src_in, K_terms+1, i,double)=x;
				MAT_RM_GET(src_in, K_terms+2, i,double)=y;
				break;
			case 3:
				MAT_RM_GET(src_in, i, K_terms+1, double)=x;
				MAT_RM_GET(src_in, i, K_terms+2, double)=y;
				MAT_RM_GET(src_in, i, K_terms+3, double)=z;
				MAT_RM_GET(src_in, K_terms+1, i, double)=x;
				MAT_RM_GET(src_in, K_terms+2, i, double)=y;
				MAT_RM_GET(src_in, K_terms+3, i, double)=z;				
				break;
		}
		

		//construct sp_src matrix(matrix that stores control points)
		switch (dim) {
			case 2:
				MAT_RM_GET(sp_src, i, 0, double)=x;
				MAT_RM_GET(sp_src, i, 1, double)=y;
				break;
			case 3:
				MAT_RM_GET(sp_src, i, 0, double)=x;
				MAT_RM_GET(sp_src, i, 1, double)=y;
				MAT_RM_GET(sp_src, i, 2, double)=z;
				break;
		}
		
	}
	
	//construct O
	for (i=0;i<dim;i++){
		for (d=0;d<dim;d++){
			MAT_RM_GET(src_in, K_terms+i, K_terms+d, double)=0.0;
		}
	}
 
	return SUCCESS;
}
 
//make the system matrix for affine
static int make_affine_matrix(Mat_rm* pts_in, Mat_rm* mat_out, const int dim) {

	int i, j;

	const int num_rows = pts_in->num_rows;

	mat_out->type = DOUBLE;
        mat_out->num_rows = num_rows;
	mat_out->num_cols = dim + 1;
        if (resize_Mat_rm(mat_out))
		return FAILURE;
        
	for (i = 0; i < num_rows; i++){

		//Add one row to the matrix
		for (j = 0; j < dim; j++) {
		    MAT_RM_GET(mat_out, i, j, double) =	
			MAT_RM_GET(pts_in, i, j, double);
		}
                MAT_RM_GET(mat_out, i, dim, double) = 1.0;
	}

	return SUCCESS;
}

//extract the control matrix from tform struct (only valid for spline)
static Mat_rm* extract_ctrl_pts(void* tform, tform_type type){
	Mat_rm* T;
	Tps *tps = (Tps* ) tform;
	switch(type){
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

static Mat_rm* extract_ctrl_pts_Tps(Tps *tps) {
	Mat_rm *kp_src = &tps->kp_src;
	return kp_src;
}

//Solve the transformation struct
static int solve_system(void* tform, Mat_rm* src, Mat_rm* ref, const int dim, 
                 tform_type type){
	//r -- the select vector to select points in source and reference matrix
	//num_pts -- number of points selected
	
	//Mat_rm *kp_ref;
	Mat_rm ref_sys, X;
	int ret;

	//extract matrices from struct
	//T = extract_params_Mat_rm(tform, type);
			
	init_Mat_rm(&ref_sys, 0, 0, DOUBLE, FALSE);
        init_Mat_rm(&X, 0, 0, DOUBLE, 0);

	//construct source matrix and initialize reference vector
	switch(type){
		case TPS: 	
                        //kp_ref = extract_ctrl_pts(tform, type);
//			make_spline_matrix(ref, &ref_in, kp_ref, num_pts, r, dim);
                        puts("solve_system: TPS not yet implemented");
                        goto SOLVE_SYSTEM_FAIL;
		case AFFINE:
			make_affine_matrix(ref, &ref_sys, dim);
			break;			
		default:
			puts("solve_system: unknown type");
			goto SOLVE_SYSTEM_FAIL;
	}	

        // solve for the coefficients					
        if (ref_sys.num_rows == ref_sys.num_cols)
          ret = solve_Mat_rm(&ref_sys, src, -1, &X);
        else 
          ret = solve_Mat_rm_ls(&ref_sys, src, &X);
        
        switch (ret) {
	    case SUCCESS:	
		break;
	    case SINGULAR:
		goto SOLVE_SYSTEM_SINGULAR;
	    default:
		goto SOLVE_SYSTEM_FAIL;
	}

	// Save the transformation matrix
	switch (type) {
		case TPS:
			//TODO
			goto SOLVE_SYSTEM_FAIL;
		case AFFINE: {

		    Mat_rm X_trans;

		    init_Mat_rm(&X_trans, 0, 0, DOUBLE, FALSE);

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

        cleanup_Mat_rm(&ref_sys);	
        cleanup_Mat_rm(&X);
	
	return SUCCESS;

SOLVE_SYSTEM_SINGULAR:	
    cleanup_Mat_rm(&ref_sys);
    cleanup_Mat_rm(&X);
    return SINGULAR;

SOLVE_SYSTEM_FAIL:
    cleanup_Mat_rm(&X);
    cleanup_Mat_rm(&ref_sys);	
    return FAILURE;
}

//Find the SSD error for the i'th point
static double tform_err_sq(void* tform, Mat_rm* src, Mat_rm* ref, int i, 
			   tform_type type){

	double err = 0.0;
	//Initialization
	//in -- inputs coordinates of source points
	//out -- registered points
	//r -- reference points (ground truth)
	double x_in, y_in, z_in, x_r, y_r, z_r, x_out, y_out, z_out;
	
	//Find the source point
	x_in = MAT_RM_GET(ref,i,0,double);
	y_in = MAT_RM_GET(ref,i,1,double);
	z_in = MAT_RM_GET(ref,i,2,double);
		
	//Register
	if(apply_tform_xyz(x_in, y_in, z_in, &x_out, &y_out, &z_out, type, 
			   tform))
		return -1.0;
		
	//Find the reference point
	x_r = MAT_RM_GET(src,i,0,double);
	y_r = MAT_RM_GET(src,i,1,double);
	z_r = MAT_RM_GET(src,i,2,double);		
	
	//Find the SSD error
	err = (x_r - x_out) * (x_r - x_out) + (y_r - y_out) * (y_r - y_out) +
		(z_r - z_out) * (z_r - z_out);	

	//return the result 
	return err;
}

//RANSAC for one iteration
static int ransac(Mat_rm* src, Mat_rm* ref, Ransac* ran, const int dim, void *tform, 
           tform_type type, int **cset, int *len) {
/* tform -- the transformation struct
	 src -- matrix of all source points [number of points * dim]
	 ref -- matrix of all reference points [number of points * dim]
*/

	/* Initialization */
	Mat_rm src_rand, ref_rand;
	int i, sel_pts, cset_len;

	const double err_thresh = ran->err_thresh;
	const double err_thresh_sq = err_thresh * err_thresh;
        const int num_src = src->num_rows;
  
	// Verify inputs
	if (src->type != DOUBLE || src->type != ref->type) {
          puts("ransac: all matrices must have type double \n");
          return FAILURE;
	}
	if (src->num_rows != ref->num_rows || src->num_cols != ref->num_cols) {
	  puts("ransac: src and ref must have the same dimensions \n");
	  return FAILURE;
	}
	
        // Initialize
        init_Mat_rm(&src_rand, 0, 0, INT, FALSE);
        init_Mat_rm(&ref_rand, 0, 0, INT, FALSE);
        
	/*Fit random points*/
	//number of points it randomly chooses
	switch(type){
		case AFFINE:
			sel_pts = dim + 1;
			break;			
		default:
			printf("ransac: unknown transformation type \n");
			goto RANSAC_FAIL;
	}

	//choose random points
	if (rand_rows(src, ref, &src_rand, &ref_rand, sel_pts))
          goto RANSAC_FAIL;
	
	//solve the system
	switch (solve_system(tform, &src_rand, &ref_rand, dim, type)) {
		case SUCCESS:
			break;
		case SINGULAR:
			goto RANSAC_SINGULAR;
		default:
			goto RANSAC_FAIL;
	}

	/*Extract consensus set*/
	//Pointwise transformation to find consensus set
	//test for each source point
	cset_len = 0;
	for (i = 0; i<num_src; i++) {

            // Calculate the error
            const double err_sq = tform_err_sq(tform, src, ref, i, type);

	    // Reject points below the error threshold
	    if (err_sq > err_thresh_sq)
		continue;

            // Add to the consensus set
            if ((*cset = realloc(*cset, ++cset_len * sizeof(int))) == NULL)
                goto RANSAC_FAIL;

            (*cset)[cset_len - 1] = i;
	}	

	// Return the new length of cset
	*len = cset_len;

        cleanup_Mat_rm(&src_rand);
        cleanup_Mat_rm(&ref_rand);
	return SUCCESS;

RANSAC_SINGULAR:
        cleanup_Mat_rm(&src_rand);
        cleanup_Mat_rm(&ref_rand);
        return SINGULAR;

RANSAC_FAIL:
        cleanup_Mat_rm(&src_rand);
        cleanup_Mat_rm(&ref_rand);
        return FAILURE;
}


//Resize spline struct based on number of selected points
int resize_Tps(Tps* tps, int num_pts, int dim){
	Mat_rm *params = &(tps->params);
	Mat_rm *kp_src = &(tps->kp_src);
	params->num_cols = num_pts + dim + 1;
	params->num_rows = dim;
	kp_src->num_rows = num_pts;
	kp_src->num_cols = dim;
	if (resize_Mat_rm(params)){
		return FAILURE;
	}
	if (resize_Mat_rm(kp_src)){
		return FAILURE;
	}
	
	tps->dim = dim;
	return SUCCESS;
}

int find_tform_ransac(Ransac* ran, Mat_rm* src, Mat_rm* ref, const int dim,
		      const tform_type type, void* tform) {

	 /* tform --transformation matrix (will fill in this struct)
		src --source points [number of points * dim]
		ref -- reference points [number of points * dim]
	 */
	
	 Mat_rm ref_cset, src_cset;	 
	 void *tform_cur;	
	 int *cset, *cset_best;
	 int i, j, num_terms, ret, len, len_best, min_num_inliers, 
	     type_min_inliers;
	 
	const int num_iter = ran->num_iter; 
        const int num_pts = src->num_rows; 	 	 	 
        const size_t tform_size = tform_get_size(type); 

	// Initialize data structures
        cset = cset_best = NULL;
        len_best = 0;

        if ((tform_cur = malloc(tform_size)) == NULL || 
	    init_tform(tform_cur, type))
          goto FIND_TFORM_FAIL;

	// initialize type-specific variables
	switch (type) {
		case AFFINE:
			num_terms = dim + 1;
			type_min_inliers = 5;
			break;
		default:
			puts("find_tform_ransac: unsupported transformation "
			     "type \n");
			return FAILURE;
	}

        min_num_inliers = (int) MAX(ceil(ran->min_inliers * num_pts), 
			      type_min_inliers);

	 if (num_pts < num_terms) {
		printf("Not enough matched points \n");
		goto FIND_TFORM_FAIL;
	 }
	
	 // Ransac iterations
	 for (i = 0; i < num_iter; i++) {
	    do {
		ret = ransac(src, ref, ran, dim, tform_cur, type, &cset, 
		             &len);
	    } while(ret == SINGULAR);

	    if (ret == FAILURE)
		goto FIND_TFORM_FAIL;

	    if (len > len_best) {
                len_best = len;
                if ((cset_best = (int *) realloc(cset_best, 
                     len * sizeof(int))) == NULL ||
                     copy_tform(tform_cur, tform, type))
		    goto FIND_TFORM_FAIL;
                memcpy(cset_best, cset, len * sizeof(int));
	    } 
	}

        // Check if the minimum number of inliers was found
        if (len_best < min_num_inliers) {
          puts("find_tform_ransac: No good model was found! \n");
          goto FIND_TFORM_FAIL;
        }

	// Initialize the concensus set matrices
	if (init_Mat_rm(&src_cset, len_best, IM_NDIMS, DOUBLE, FALSE) ||
	    init_Mat_rm(&ref_cset, len_best, IM_NDIMS, DOUBLE, FALSE))
	    goto FIND_TFORM_FAIL;

	// extract the concensus set
	MAT_RM_LOOP_START(&src_cset, i, j)
	    
	    const int idx = cset_best[i];

	    MAT_RM_GET(&src_cset, i, j, double) = 
		MAT_RM_GET(src, idx, j, double);
	    MAT_RM_GET(&ref_cset, i, j, double) = 
		MAT_RM_GET(ref, idx, j, double);

	MAT_RM_LOOP_END

#ifdef RANSAC_REFINE
	// Refine with least squares
	switch (solve_system(tform_cur, &src_cset, &ref_cset, dim, type)) {
	    case SUCCESS:
		// Copy the refined transformation to the output
		if (copy_tform(tform_cur, tform, type))
		    goto FIND_TFORM_FAIL;
		break;
	    case SINGULAR:
		// Stick with the old transformation 
#ifdef VERBOSE
		printf("find_tform_ransac: warning: least-squares refinement "
		       "abandoned due to numerical precision \n");
#endif
		break;
	   default:
		goto FIND_TFORM_FAIL; 
	}
#endif
	 
        free(cset);
        free(cset_best);
        cleanup_tform(tform_cur, type);
	if (tform_cur != NULL)
	    free(tform_cur);
	return SUCCESS;

FIND_TFORM_FAIL:
        if (cset != NULL)
          free(cset);
        if (cset_best != NULL)
          free(cset_best);
        cleanup_tform(tform_cur, type);
	if (tform_cur != NULL)
	    free(tform_cur);
	return FAILURE;
 }

 
 
