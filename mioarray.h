#ifndef ARRAY_H
#define ARRAY_H

#include "mio.h"

struct array_stc_2d
{
        int xdim;
        int ydim;
	MIO *mio;
        double *data;
};

typedef struct array_stc_2d Array2D;
Array2D *MakeArray2D(int xdim, int ydim);
void FreeArray2D(Array2D *a);
void PrintArray2D(Array2D *a);
Array2D *CopyArray2D(Array2D *a);
Array2D *TransposeArray2D(Array2D *a);
Array2D *MultiplyArray2D(Array2D *a, Array2D *b);
Array2D *InvertArray2D(Array2D *a);

/*
 * choose 1D to be column vectors
 */

#define Array1D Array2D
#define MakeArray1D(dim) MakeArray2D((dim),1)
#define FreeArray1D(a) FreeArray2D(a)
#define PrintArray1D(a) PrintArray2D(a)
#define CopyArray1D(a) CopyArray2D(a)
#define TransposeArray1D(a) TransposeArray2D(a)

#endif

