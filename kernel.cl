#pragma OPENCL EXTENSION cl_khr_fp64: enable
#include "kernel.h"

__kernel void clkernel(const size_t numGroups,
                       global scoreType *scores,
                       global  blockOffsetType * blockOffsets,
                        global seqNumType* seqNums,
                        global seqType* sequences
                      //  global unsigned short *matrix_in,
                        //global unsigned short *vector_out
){
    scoreType st;
    int groupNum = get_global_id(0);
  //  scores[xx] =  20;
    if (groupNum == 2){ // 2
        int seqBlock = groupNum >> 4;  //GPUdb::LOG2_BLOCK_SIZE;
        int groupNumInBlock = groupNum % 15; //BLOCK_SIZE - 1
        int groupOffset = blockOffsets[seqBlock]+(groupNumInBlock*8); //SUBBLOCK_SIZE = 8
        //printf("d %d\n", seqBlock);
       // seqType* sequence = sequences[0];
        __global seqType* sequence;
        sequence = &sequences[groupOffset];
        seqType8 s = *(seqType8*)sequence;
        int i = 0;
       while(s.a.x!=' ') //Until terminating subblock
        {
            if(s.a.x=='#') //Subblock signifying concatenated sequences
            {
                printf("\n");
            }
    //        printf("%d %d %d %d ", s.a.x, s.a.y, s.a.z, s.a.w);
      //      printf("%d %d %d %d ", s.b.x, s.b.y, s.b.z, s.b.w);
           // printf("%d %d %d %d\n", s.b.w, s.b.x, s.b.y, s.b.z);
            sequence += 16*8;
            s = *(seqType8*)sequence;
        }
        printf("\n---\n");
}
}