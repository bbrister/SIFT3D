/* diffeoSift3d.c 
 * ----------------------------------------------------------------
 * Extract dense SIFT3D features from an image, and perform
 * diffeomorphic registration.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "macros.h"
#include "imutil.h"
#include "sift.h"
#include "reg.h"

#define REG_PATH "out/reg.nii"
#define GRID_PATH "out/grid.nii"
#define SRC_OUT_PATH "out/src.nii"
#define REF_OUT_PATH "out/ref.nii"

int main(int argc, char **argv) {

        Image src, ref, srcp, refp, reg, desc_src, desc_ref, tform, grid, grid_reg;
        SIFT3D sift3d;
        char *src_path, *ref_path;
        int nx, ny, nz;

        /* Parse the arguments */
        if (argc != 3) {
                puts("Usage: denseSift3d [source.nii] [reference.nii] \n");
                return 1;
        }
        src_path = argv[1];
        ref_path = argv[2];

        /* Initialize data */
        init_im(&src);
        init_im(&ref);
        init_im(&reg);
        init_im(&desc_src);
        init_im(&desc_ref);
        init_im(&tform);
        if (init_SIFT3D(&sift3d))
                err_exit("initalize sift3d \n");

        /* Read the images */        
        if (read_nii(src_path, &src))
                err_exit("read the source image \n");                
        if (read_nii(ref_path, &ref))
                err_exit("read the source image \n");                

        /* Pad the images to have the same size */
	nx = MAX(src.nx, ref.nx);
	ny = MAX(src.ny, ref.ny);
	nz = MAX(src.nz, ref.nz);
	if (init_im_first_time(&srcp, nx, ny, nz, 1) || 
	    init_im_first_time(&refp, nx, ny, nz, 1) || 
	    im_pad(&src, &srcp) || 
	    im_pad(&ref, &refp))
	    err_exit("pad images \n");

        /* Release the un-padded images */
        im_free(&src);
        im_free(&ref);

        /* Write the inputs */
        if (write_nii(SRC_OUT_PATH, &srcp))
                err_exit("write the source image \n");                
        if (write_nii(REF_OUT_PATH, &refp))
                err_exit("write the source image \n");                
#if 1
        /* Extract the descriptors */
        if (SIFT3D_extract_dense_descriptors(&sift3d, &srcp, &desc_src))
                err_exit("extract source descriptors \n");
        if (SIFT3D_extract_dense_descriptors(&sift3d, &refp, &desc_ref))
                err_exit("extract source descriptors \n");
#else
        desc_src = srcp;
        desc_ref = refp;
#endif

        /* Register */
        if (diffeo_reg(&desc_src, &desc_ref, &tform))
                err_exit("register \n");

        /* Transform the source image */
        if (im_inv_transform(&srcp, &reg, WARP, &tform, LINEAR))
                err_exit("transform source \n");

        /* Draw a grid */
        init_im(&grid);
	if (draw_grid(&grid, srcp.nx, srcp.ny, srcp.nz, 20, 1))
		err_exit("make grid \n");

        /* Transform the grid */
        init_im(&grid_reg);
        if (im_inv_transform(&grid, &grid_reg, WARP, &tform, LINEAR))
                err_exit("transform grid \n");

        /* Write the output images */
        if (write_nii(REG_PATH, &reg))
                err_exit("write reg\n");
        if (write_nii(GRID_PATH, &grid_reg))
                err_exit("write grid\n");

        return 0;
}
