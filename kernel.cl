//#pragma OPENCL EXTENSION cl_khr_fp64: enable
#include "kernel.h"

__kernel void clkernel(const size_t numGroups,
                       global scoreType *scores,
                       global  blockOffsetType * blockOffsets
                      //  global unsigned short *matrix_in,
                        //global unsigned short *vector_out
){
    scoreType st;

    int xx = get_global_id(0);
    scores[xx] =  20;
    if (xx == 0){
        //printf("%d\n", numGroups);
        printf("d %d\n", blockOffsets[1]);
        printf("d %d\n", blockOffsets[2]);
        printf("d %d\n", blockOffsets[3]);
    }
 //   printf("%d\n", numGroups);
 //vector_out[x] = 20;
}