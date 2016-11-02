//
// Created by Bjorn Sigurbergsson on 01/11/16.
//

typedef unsigned short scoreType;
typedef unsigned int blockOffsetType;
typedef unsigned int seqSizeType;
typedef unsigned int seqNumType;

struct TempData
{
    scoreType F;
    scoreType Ix;
};

struct TempData2
{
    struct TempData a;
    struct TempData b;
};