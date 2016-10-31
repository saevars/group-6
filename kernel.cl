#pragma OPENCL EXTENSION cl_khr_fp64: enable
__kernel void clkernel(global char *vector_in, global unsigned short *matrix_in, global unsigned short *vector_out){
 int x = get_global_id(0);
 vector_out[x] = 20;
}