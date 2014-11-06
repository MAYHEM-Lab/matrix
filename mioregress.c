#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mioarray.h"


/*
 * normal equations are
 * X^tXB = X^ty
 *
 * B = (X^tX)^-1 * X^t * y
 */
Array1D *Regress(Array2D *x, Array2D *y)
{
	Array1D *B;
	Array2D *xt; /* x^t */
	Array2D *xtx; /* x^t * x */
	Array2D *xtxI; /* (x^t * x)^-1 */
	Array2D *xtxIxt; /* ((x^t * x)^-1) * x^t */
	int i;

	xt = TransposeArray2D(x);
	if(xt == NULL)
	{
		return(NULL);
	}

	xtx = MultiplyArray2D(xt,x);
	if(xtx == NULL)
	{
		FreeArray2D(xt);
		return(NULL);
	}

	xtxI = InvertArray2D(xtx);
	FreeArray2D(xtx);
	if(xtxI == NULL)
	{
		FreeArray2D(xt);
		return(NULL);
	}

	xtxIxt = MultiplyArray2D(xtxI,xt);
	FreeArray2D(xt);
	FreeArray2D(xtxI);
	if(xtxIxt == NULL)
	{
		return(NULL);
	}
	B = MultiplyArray2D(xtxIxt,y);
	FreeArray2D(xtxIxt);
	if(B == NULL)
	{
		return(NULL);
	}

	return(B);
}
	

#ifdef STANDALONE

int main(int argc, char *argv[])
{
	int i;
	int j;
	Array2D *a;
	Array2D *b;
	Array2D *c;

	a = MakeArray2D(7,2);
	b = MakeArray2D(7,1);

	a->data[0*2+0] = 1;
	a->data[0*2+1] = 1232;
	a->data[1*2+0] = 1;
	a->data[1*2+1] = 1070;
	a->data[2*2+0] = 1;
	a->data[2*2+1] = 1086;
	a->data[3*2+0] = 1;
	a->data[3*2+1] = 1287;
	a->data[4*2+0] = 1;
	a->data[4*2+1] = 1130;
	a->data[5*2+0] = 1;
	a->data[5*2+1] = 1048;
	a->data[6*2+0] = 1;
	a->data[6*2+1] = 1121;

	b->data[0] = 3.52;
	b->data[1] = 2.91;
	b->data[2] = 2.4;
	b->data[3] = 3.47;
	b->data[4] = 3.47;
	b->data[5] = 2.37;
	b->data[6] = 2.4;


	c = Regress(a,b);
	if(c == NULL) {
		fprintf(stderr,"regress 1 failed\n");
		exit(1);
	}

	printf("a: ");
	PrintArray2D(a);
	printf("b: ");
	PrintArray2D(b);
	printf("c: ");
	PrintArray2D(c);

	FreeArray2D(a);
	FreeArray2D(b);
	FreeArray2D(c);
	return(0);
}

#endif

	

	



