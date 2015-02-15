/* reg.h
 * ----------------------------------------------------------------
 * Rice MRI Team
 * ----------------------------------------------------------------
 * Internal header for reg.c
 *-----------------------------------------------------------------
 * Created: Jenny Zhang 1/24/2014
 * Last updated: Blaine Rister 11/18/2014
 */
#ifndef _REG_H
#define _REG_H

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "types.h"
#include "macros.h"
#include "imutil.h"

/* Implementaiton parameters */
#define RANSAC_REFINE	// Use least-squares refinement in RANSAC

static const double singular_thresh = 1e-15; 
static const double err_init = 1e40; 
static const double spline_percent = 0.45; 

int init_Ransac(Ransac *ran, double min_inliers, double err_thresh, 
				 int num_iter);
					  
int find_tform_ransac(Ransac* ran, Mat_rm* src, Mat_rm* ref, const int dim,
					  tform_type type, void* tform);
					  					  
#endif
