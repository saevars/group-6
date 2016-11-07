#pragma OPENCL EXTENSION cl_khr_fp64: enable
#include "kernel.h"

constant short BLOCK_SIZE = 16;
constant short SUBBLOCK_SIZE = 8;
constant short LOG2_BLOCKSIZE = 4;
constant int queryProfileLength = 320;
constant int noOfThreads = 1024;
constant int gapPenaltyTotal = -12;
constant int gapExtendPenalty = -2;

char4 populateSubstScoreFromQueryProfile(global substType *queryProfile, char a, int j){
    char4 score;
    score.x = queryProfile[a * queryProfileLength + 4*j];
    score.y = queryProfile[a * queryProfileLength + 4*j+1];
    score.z = queryProfile[a * queryProfileLength + 4*j+2];
    score.w = queryProfile[a * queryProfileLength + 4*j+3];
    return score;
}

/*void alignResidues(scoreType maxScore, scoreType4 left, int4 ixLeft, scoreType top,
                   scoreType topLeft, int IyTop, char4 substScore)
{

}*/

residue alignResidues(residue res, char4 substScore){//, char4 substScores){
    //q0
    res.ixLeft.x = max(0,max(res.left.x+gapPenaltyTotal,res.ixLeft.x+gapExtendPenalty)); //Max(0,...) here so IxColumn[] can be unsigned
    res.IyTop = max(res.top+gapPenaltyTotal,res.IyTop+gapExtendPenalty);
     int align = res.topLeft+substScore.x;
     res.topLeft=res.left.x;
     res.left.x = max(align,max(res.ixLeft.x,res.IyTop));
    //q1
     res.ixLeft.y = max(0,max(res.left.y+gapPenaltyTotal,res.ixLeft.y+gapExtendPenalty)); //Max(0,...) here so IxColumn[] can be unsigned
     res.IyTop = max(res.left.x+gapPenaltyTotal,res.IyTop+gapExtendPenalty);
     align = res.topLeft+substScore.y;
     res.topLeft=res.left.y;
     res.left.y = max(align,max(res.ixLeft.y,res.IyTop));

     //q2
     res.ixLeft.z = max(0,max(res.left.z+gapPenaltyTotal,res.ixLeft.z+gapExtendPenalty)); //Max(0,...) here so IxColumn[] can be unsigned
     res.IyTop = max(res.left.y+gapPenaltyTotal,res.IyTop+gapExtendPenalty);
     align = res.topLeft+substScore.z;
     res.topLeft=res.left.z;
     res.left.z = max(align,max(res.ixLeft.z,res.IyTop));

     //q3
     res.ixLeft.w = max(0,max(res.left.w+gapPenaltyTotal,res.ixLeft.w+gapExtendPenalty)); //Max(0,...) here so IxColumn[] can be unsigned
     res.IyTop = max(res.left.z+gapPenaltyTotal,res.IyTop+gapExtendPenalty);
     align = res.topLeft+substScore.w;
     res.left.w = max(align,max(res.ixLeft.w,res.IyTop));

     res.topLeft=res.top; //The next column is to the right of this one, so current res.top res.left becomes new res.top
     res.top = res.left.w; //Set res.top value for next query chunk
     res.maxScore = max(res.left.x,max(res.left.y,max(res.left.z,max(res.left.w,res.maxScore)))); //Update max score
     return res;
}

/*
 * alignWithQuery((substType *)queryProfile, (seqType8) s,
                       (TempData2*) tempColumn, (scoreType) maxScore, column);*/
scoreType alignWithQuery(global substType *queryProfile, char8 s,
                    global TempData2*  tempColumn, scoreType maxScore, int column, seqNumType seqNum){

    //Set the top related values to 0 as we're at the top of the matrix
    scoreType8 top = {0,0,0,0,0,0,0,0};
    scoreType topLeft = 0;
    int8 IyTop = {0,0,0,0,0,0,0,0};

    char4 substScores; //Query profile scores
    scoreType4 left;
    int4 ixLeft;

    residue res, res2;


    for(int j = 0; j < queryProfileLength/4; j++)
    {
        if(seqNum == 16){
      //      printf("j is %d and tempColumn %d %d \n", j, tempColumn[0].a, tempColumn[0].b);
        }
        TempData2 t = tempColumn[0];
        left.x = column * t.a.F;
        ixLeft.x = column* t.a.Ix;
        left.y = column*t.b.F;
        ixLeft.y = column* t.b.Ix;

        //Load second half of temporary column
        t = tempColumn[noOfThreads];
        left.z = column*t.a.F;
        ixLeft.z = column*t.a.Ix;
        left.w = column*t.b.F;
        ixLeft.w = column*t.b.Ix;

        int topLeftNext = left.w; //Save the top left cell value for the next loop interation


        res.maxScore = maxScore;


        res.left = left;
        res.ixLeft = ixLeft;
        res.topLeft = topLeft;


        res.top = (scoreType)top.a.x;
        res.IyTop = IyTop.lo.x;

//        if (j == 0 && seqNum == 15){
//            printf("j = 1; after maxscore = %d\n",
//                   res.maxScore);
//        }

        substScores = (char4) populateSubstScoreFromQueryProfile(queryProfile, s.lo.x, j);
        res2 = alignResidues((residue)res, (char4) substScores);

        top.a.x = (scoreType)res2.top;
        IyTop.lo.x = res2.IyTop;

        res2.top = (scoreType)top.a.y;
        res2.IyTop = IyTop.lo.y;

        substScores = (char4) populateSubstScoreFromQueryProfile(queryProfile, s.lo.y, j);
        res = alignResidues((residue)res2, (char4) substScores);

        top.a.y = (scoreType)res.top;
        IyTop.lo.y = res.IyTop;

        res.top = (scoreType)top.a.z;
        res.IyTop = IyTop.lo.z;

        substScores = (char4) populateSubstScoreFromQueryProfile(queryProfile, s.lo.z, j);
        res2 = alignResidues((residue)res, (char4) substScores);

        top.a.z = (scoreType)res2.top;
        IyTop.lo.z= res2.IyTop;

        res2.top = (scoreType)top.a.w;
        res2.IyTop = IyTop.lo.w;

        substScores = (char4) populateSubstScoreFromQueryProfile(queryProfile, s.lo.w, j);
        res = alignResidues((residue)res2, (char4) substScores);

        top.a.w = (scoreType)res.top;
        IyTop.lo.w = res.IyTop;


        res.top = (scoreType)top.b.x;
        res.IyTop = IyTop.hi.x;

        substScores = (char4) populateSubstScoreFromQueryProfile(queryProfile, s.hi.x, j);
        res2 = alignResidues((residue)res, (char4) substScores);

        top.b.x = (scoreType)res2.top;
        IyTop.hi.x = res2.IyTop;



        res2.top = (scoreType)top.b.y;
        res2.IyTop = IyTop.hi.y;

        substScores = (char4) populateSubstScoreFromQueryProfile(queryProfile, s.hi.y, j);
        res = alignResidues((residue)res2, (char4) substScores);

        top.b.y = (scoreType)res.top;
        IyTop.hi.y = res.IyTop;




        res.top = (scoreType)top.b.z;
        res.IyTop = IyTop.hi.z;

        substScores = (char4) populateSubstScoreFromQueryProfile(queryProfile, s.hi.z, j);
        res2 = alignResidues((residue)res, (char4) substScores);

        top.b.z = (scoreType)res2.top;
        IyTop.hi.z = res2.IyTop;
        /*
         * ********
         */

        res2.top = (scoreType)top.b.w;
        res2.IyTop = IyTop.hi.w;

        /*if (j == 1){
            printf("before maxscore = %d, left = {%d, %d, %d, %d}, ixLeft = {%d, %d, %d, %d}, top = %d, topLeft = %d, IyTop = %d\n",
                   res2.maxScore, res2.left.x, res2.left.y, res2.left.z, res2.left.w,
                   res2.ixLeft.x,res2.ixLeft.y,res2.ixLeft.z,res2.ixLeft.w,
                   res2.top, res2.topLeft, res2.IyTop);
        }*/

        substScores = (char4) populateSubstScoreFromQueryProfile(queryProfile, s.hi.w, j);
        res = alignResidues((residue)res2, (char4) substScores);

        top.b.w = (scoreType)res.top;
        IyTop.hi.w = res.IyTop;




        topLeft = topLeftNext;
        maxScore = res.maxScore;

        //Save the two temporary column values
        t.a.F = res.left.x;

        t.a.Ix = res.ixLeft.x;
        t.b.F = res.left.y;
        t.b.Ix = res.ixLeft.y;
        tempColumn[0]=t;
        tempColumn += noOfThreads;
        t.a.F = res.left.z;
        t.a.Ix = res.ixLeft.z;
        t.b.F = res.left.w;
        t.b.Ix = res.ixLeft.w;
        tempColumn[0]=t;
     //   if (j == 0)printf("%d %d\n", tempColumn[0].a.Ix, t.a.Ix);

        tempColumn += noOfThreads;

        if (j < 3){
          // printf("%d\n", (int)res.maxScore);
         //   printf("%d, %d, %d, %d\n", substScores.x, substScores.y, substScores.z, substScores.w);
        }
    }
    return maxScore;
}


void align(global seqType* sequence, global TempData2* tempColumn, seqNumType seqNum,
          global scoreType* scores, global substType *queryProfile){

    scoreType maxScore=0;
    int column = 0; //Column = 0 means that alignment function will use 0 for 'left' values as there's no left column to read from

    char8 s ;
    __global seqType * tempSeq = sequence;
    s.lo.x = * tempSeq++;
    s.lo.y = * tempSeq++;
    s.lo.z = * tempSeq++;
    s.lo.w = * tempSeq++;
    s.hi.x = * tempSeq++;
    s.hi.y = * tempSeq++;
    s.hi.z = * tempSeq++;
    s.hi.w = * tempSeq++;

    if (seqNum == 0){
        printf("s.lo.x = %d \n", s.lo.x);
        printf("s.lo.y = %d \n", s.lo.y);
        printf("s.lo.z = %d \n", s.lo.z);
        printf("s.lo.w = %d \n", s.lo.w);
        printf("s.hi.x = %d \n", s.hi.x);
        printf("s.hi.y = %d \n", s.hi.y);
        printf("s.hi.z = %d \n", s.hi.z);
        printf("s.hi.w = %d \n\n", s.hi.w);
}

    while(s.lo.x!=' ') //Until terminating subblock
    {
        if(s.lo.x=='#') //Subblock signifying concatenated sequences
        {
            scores[seqNum] = maxScore; //Set score for sequence
            seqNum++;
            column=maxScore=0;
        }
        maxScore = alignWithQuery(queryProfile, s,
                        tempColumn, (scoreType) maxScore, column, (seqNumType) seqNum);

        column=1;
        sequence += BLOCK_SIZE*SUBBLOCK_SIZE;

        tempSeq = sequence;
        s.lo.x = * tempSeq++;
        s.lo.y = * tempSeq++;
        s.lo.z = * tempSeq++;
        s.lo.w = * tempSeq++;
        s.hi.x = * tempSeq++;
        s.hi.y = * tempSeq++;
        s.hi.z = * tempSeq++;
        s.hi.w = * tempSeq++;
    }

        scores[seqNum] = maxScore;
}


/*printf("%d %d %d %d ", s.a.x, s.a.y, s.a.z, s.a.w);
        printf("%d %d %d %d ", s.b.x, s.b.y, s.b.z, s.b.w);
        printf("%d %d %d %d\n", s.b.w, s.b.x, s.b.y, s.b.z);*/
__kernel void clkernel(const unsigned int numGroups,
                      global scoreType *scores,
                      global  blockOffsetType * blockOffsets,
                      global  seqNumType* seqNums,
                      global  seqType* sequences,
                      global  substType *queryProfile,
                       global TempData2 *tempColumns
){

    scoreType st;
    int groupNum = get_global_id(0);
    __global TempData2* tempColumn = &tempColumns[groupNum];

  //  scores[xx] =  20;
    if (groupNum < numGroups){
        seqNumType seqNum=seqNums[groupNum];

        int seqBlock = groupNum >> LOG2_BLOCKSIZE;
        int groupNumInBlock = groupNum % (BLOCK_SIZE);
        int groupOffset = blockOffsets[seqBlock]+(groupNumInBlock*SUBBLOCK_SIZE);

       // printf("seq=%d, seqBlock=%d, gnib=%d, goffset=%d\n",
         //   groupNum, seqBlock, groupNumInBlock, groupOffset);


__global seqType* sequence = &sequences[groupOffset];

        align(
         sequence, tempColumn, seqNum, scores,
         queryProfile);
        groupNum += noOfThreads;
    /*
     * void align(seqType* sequence, TempData2* const tempColumn,
           seqSizeType& seqNum, scoreType* scores)*/

}
}

