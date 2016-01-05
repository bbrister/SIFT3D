/* -----------------------------------------------------------------------------
 * kernels.cl 
 * -----------------------------------------------------------------------------
 * Copyright (c) 2015-2016 Blaine Rister et al., see LICENSE for details.
 * -----------------------------------------------------------------------------
 * This file contains OpenCL kernels for image processing.
 * -----------------------------------------------------------------------------
 */

#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

sampler_t sampler_downsample_2x = CLK_NORMALIZED_COORDS_FALSE |
				  CLK_ADDRESS_CLAMP_TO_EDGE |
				  CLK_FILTER_NEAEREST;

kernel void downsample_2x_3d(__read_only image3d_t src,
						 	 __write_only image3d_t dst) {

	int x, y, z;
	float4 out;

	x = get_global_id(0);
	y = get_global_id(1);
	z = get_global_id(2);
	out = read_imagef(src, sampler_downsample_2x, (int4) (x, y, z, 0) * 2);
	write_imagef(dst, (int4) (x, y, z, 0), out);
}
