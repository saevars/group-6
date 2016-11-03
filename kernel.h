//
// Created by Bjorn Sigurbergsson on 01/11/16.
//

typedef unsigned short scoreType;
typedef unsigned int blockOffsetType;
typedef signed char seqType;
typedef unsigned int seqNumType;
//typedef char4 seqType4;
//typedef char8 seqType8;

//typedef cl_char8 seqType8;

typedef char4 seqType4;

typedef struct seqType8
{
    seqType4 a;
    seqType4 b;

} seqType8;


/*TBA
typedef struct TempData
{
    scoreType F;
    scoreType Ix;
} TempData;

typedef struct TempData2
{
    struct TempData a;
    struct TempData b;
} TempData2;
*/