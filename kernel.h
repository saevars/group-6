//
// Created by Bjorn Sigurbergsson on 01/11/16.
//
//
typedef unsigned short scoreType;
typedef unsigned int blockOffsetType;
typedef signed char seqType;
typedef unsigned int seqNumType;
typedef ushort4 scoreType4;
////typedef char4 seqType4;size_t
////typedef char8 seqType8;
//
////typedef cl_char8 seqType8;
//
typedef char4 seqType4;
typedef char substType;

/*alignResidues(scoreType &maxScore, scoreType4& left, int4& ixLeft, scoreType& top,
	scoreType& topLeft, int &IyTop, const char4 &substScore)*/
//
typedef struct residue {
    scoreType maxScore;
    scoreType4 left;
    int4 ixLeft;
    scoreType top;
    scoreType topLeft;
    int IyTop;
} residue;

typedef struct seqType8
{
    seqType4 a;
    seqType4 b;

} seqType8;

typedef struct  scoreType8
{
    scoreType4 a;
    scoreType4 b;
} scoreType8;

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
