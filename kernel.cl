#pragma OPENCL EXTENSION cl_khr_fp64: enable
#include "kernel.h"

constant short BLOCK_SIZE = 16;
constant short SUBBLOCK_SIZE = 8;
constant short LOG2_BLOCKSIZE = 4;
constant int queryProfileLength = 256;

char4 populateSubstScoreFromQueryProfile(substType *queryProfile, char a, int j){
    char4 score;
    score.x = queryProfile[a * queryProfileLength + 4*j];
    score.y = queryProfile[a * queryProfileLength + 4*j+1];
    score.z = queryProfile[a * queryProfileLength + 4*j+2];
    score.w = queryProfile[a * queryProfileLength + 4*j+3];
    return score;
}

void alignWithQuery(substType *queryProfile, seqType8 s){
    char4 substScores;
    for(int j = 0; j < queryProfileLength; j++)
    {
        substScores = populateSubstScoreFromQueryProfile(queryProfile, s.a.x, j);

        if (j < 3){
           printf("j=%d s.a.x = %d\n", j, s.a.x);
            printf("%d, %d, %d, %d\n", substScores.x, substScores.y, substScores.z, substScores.w);
        }
    }

}

__kernel void clkernel(const size_t numGroups,
                       global scoreType *scores,
                       global  blockOffsetType * blockOffsets,
                        global seqNumType* seqNums,
                        global seqType* sequences,
                        global substType *queryProfile
                        //global unsigned short *vector_out
){
    scoreType st;
    int groupNum = get_global_id(0);
  //  scores[xx] =  20;
    if (groupNum == 2){ // 2
        int seqBlock = groupNum >> LOG2_BLOCKSIZE;
        int groupNumInBlock = groupNum % (BLOCK_SIZE - 1); //BLOCK_SIZE - 1
        int groupOffset = blockOffsets[seqBlock]+(groupNumInBlock*SUBBLOCK_SIZE); //SUBBLOCK_SIZE = 8
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
        alignWithQuery((substType *)queryProfile, (seqType8) s);
            /*printf("%d %d %d %d ", s.a.x, s.a.y, s.a.z, s.a.w);
            printf("%d %d %d %d ", s.b.x, s.b.y, s.b.z, s.b.w);
            printf("%d %d %d %d\n", s.b.w, s.b.x, s.b.y, s.b.z);*/
            sequence += BLOCK_SIZE*SUBBLOCK_SIZE;
            s = *(seqType8*)sequence;
        }
}
}