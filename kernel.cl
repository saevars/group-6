#pragma OPENCL EXTENSION cl_khr_fp64: enable
#include "kernel.h"

__kernel void clkernel(global scoreType *vector_in, global unsigned short *matrix_in, global unsigned short *vector_out){
    scoreType st;
 //int x = get_global_id(0);
 //vector_out[x] = 20;
}