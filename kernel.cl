//#pragma OPENCL EXTENSION cl_khr_fp64: enable
#include "kernel.h"

__kernel void clkernel(const size_t numGroups, global scoreType *scores//,
                      //  global unsigned short *matrix_in,
                        //global unsigned short *vector_out
){
    scoreType st;

    int xx = get_global_id(0);
   // scores[xx] = (scoreType) 20;
    printf("%d\n", numGroups);
 //vector_out[x] = 20;
}