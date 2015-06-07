#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mioarray.h"


double RSquared(Array2D *x, Array2D *b, Array2D *y)
{
	Array2D *f;
	int i;
	double ss_res;
	double ss_reg;
	double ss_total;
	double y_bar;
	double temp;
	double rsq;

	f = MultiplyArray2D(x,b);
	if(f == NULL) {
		return(-1.0);
	}

	ss_res = 0;
	for(i=0; i < f->ydim; i++) {
		temp = (y->data[i*y->xdim+0] -
			   f->data[i*f->xdim+0]);
		ss_res += (temp*temp);
	}

	y_bar = 0;
	for(i=0; i < y->ydim; i++) {
		y_bar += y->data[i*y->xdim+0];
	}
	y_bar = y_bar / (double)(y->ydim);

	ss_reg = 0;
	for(i=0; i < f->ydim; i++) {
		temp = (y_bar -
			   f->data[i*f->xdim+0]);
		ss_reg += (temp*temp);
	}

	ss_total = ss_reg + ss_res;

	rsq = 1 - (ss_res / ss_total);

	FreeArray2D(f);

	return(rsq);
}

double RMSE(Array2D *x, Array2D *b, Array2D *y)
{
	Array2D *f;
	int i;
	double ss_res;
	double temp;
	double rmse;
	double count;

	f = MultiplyArray2D(x,b);
	if(f == NULL) {
		return(-1.0);
	}

	ss_res = 0;
	count = 0;
	for(i=0; i < f->ydim; i++) {
		temp = (y->data[i*y->xdim+0] -
			   f->data[i*f->xdim+0]);
		ss_res += (temp*temp);
		count++;
	}

	rmse = sqrt(ss_res/count);

	FreeArray2D(f);

	return(rmse);
}

Array2D *Resduals(Array2D *x, Array2D *b, Array2D *y)
{
	Array2D *f;
	int i;
	double temp;
	Array2D *resid;

	f = MultiplyArray2D(x,b);
	if(f == NULL) {
		return(NULL);
	}

	resid = MakeArray1D(f->ydim);
	if(resid == NULL) {
		FreeArray2D(f);
		return(NULL);
	}

	for(i=0; i < f->ydim; i++) {
		temp = (y->data[i*y->xdim+0] -
			   f->data[i*f->xdim+0]);
		resid->data[i*resid->xdim+0] = temp;
	}


	FreeArray2D(f);

	return(resid);
}

double RSquaredOld(Array2D *x, Array2D *b, Array2D *y)
{
	Array2D *f;
	int i;
	double ss_res;
	double ss_total;
	double y_bar;
	double temp;
	double rsq;

	f = MultiplyArray2D(x,b);
	if(f == NULL) {
		return(-1.0);
	}

	ss_res = 0;
	for(i=0; i < f->ydim; i++) {
		temp = (y->data[i*y->xdim+0] -
			   f->data[i*f->xdim+0]);
		ss_res += (temp*temp);
	}

	y_bar = 0;
	for(i=0; i < y->ydim; i++) {
		y_bar += y->data[i*y->xdim+0];
	}
	y_bar = y_bar / (double)(y->ydim);

	ss_total = 0;
	for(i=0; i < y->ydim; i++) {
		temp = (y->data[i*y->xdim+0] -
			y_bar);
		ss_total += (temp*temp);
	}

	rsq = 1 - (ss_res / ss_total);

	FreeArray2D(f);

	return(rsq);
}

Array2D *CenterScale(Array2D *x)
{
	Array2D *sx;
	double mu;
	double sigma;
	double acc;
	double count;
	int i;
	int j;

	sx = MakeArray2D(2,x->xdim);
	if(sx == NULL) {
		return(NULL);
	}

	for(j=0; j < x->xdim; j++) {
		acc = 0;
		count = 0;
		for(i=0; i < x->ydim; i++) {
			acc += x->data[i*x->xdim+j];
			count++;
		}
		mu = acc / count;
		acc = 0;
		count = 0;
		for(i=0; i < x->ydim; i++) {
			acc += ((x->data[i*x->xdim+j] - mu) *
				(x->data[i*x->xdim+j] - mu));
			count++;
		}
		sigma = sqrt(acc / count);
		/*
		 * first row is mu
		 * second row is sigma
		 */
		sx->data[0*sx->xdim+j] = mu;
		sx->data[1*sx->xdim+j] = sigma;
	}

	return(sx);
}



/*
 * normal equations are
 * X^tXB = X^ty
 *
 * B = (X^tX)^-1 * X^t * y
 */
Array1D *RegressMatrix2DSimple(Array2D *x, Array2D *y)
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

#ifndef USELAPACK

Array1D *RegressMatrix2D(Array2D *x, Array2D *y)
{
	return(RegressMatrix2DSimple(x,y));
}

#else

#include "lapacke.h"

/*
 * from http://www.netlib.org/lapack/lapacke.html#_calling_code_dgels_code
 */

Array1D *RegressMatrix2D(Array2D *x, Array2D *y)
{
	Array1D *B;
	int i;
	char trans = 'N';
	lapack_int m;
	lapack_int n;
	lapack_int nrhs;
	lapack_int lda;
	lapack_int ldb;
	lapack_int info;
	int ndx;

	m = x->ydim;
	n = x->xdim;
	nrhs = 1;

	/*
	 * in row major order, lda is the number of values that separate two elements in the same
	 * column
	 */
	lda = n;
	ldb = y->xdim;

	B = MakeArray1D(y->ydim);
	if(B == NULL) {
		return(NULL);
	}

	/*
	 * dgels destroys the input B vector
	 */
	for(i=0; i < B->ydim; i++) {
		B->data[i] = y->data[i];
	}

	info = LAPACKE_dgels(LAPACK_ROW_MAJOR,'N',m,n,nrhs,
					x->data,lda,B->data,ldb);
	if(info != 0) {
		fprintf(stderr,"LAPACK error in regression: %d\n",
			info);
		ndx = fabs((double)info) - 1;
		fprintf(stderr,"x[%d]: %f, y[%d]: %f\n",
			ndx,x->data[ndx],ndx,y->data[ndx]);
		FreeArray1D(B);
		return(NULL);
	}

	/*
            On exit, if INFO = 0, B is overwritten by the solution   
            vectors, stored columnwise:   
            if TRANS = 'N' and m >= n, rows 1 to n of B contain the least   
            squares solution vectors; the residual sum of squares for the   
            solution in each column is given by the sum of squares of   
            elements N+1 to M in that column;
	*/

	B->ydim = n;

	return(B);
}
#endif
	

#ifdef STANDALONE

int main(int argc, char *argv[])
{
	int i;
	int j;
	Array2D *a;
	Array2D *ac;
	Array2D *b;
	Array2D *c;
	Array2D *d;

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


	c = RegressMatrix2D(a,b);
	if(c == NULL) {
		fprintf(stderr,"regress 1 failed\n");
		exit(1);
	}

	/*
	 * answer should be
	 * c: 	[-2.417108]
	 *      [0.004698]
	 */


	printf("a: ");
	PrintArray2D(a);
	printf("b: ");
	PrintArray2D(b);
	printf("c: ");
	PrintArray2D(c);

	printf("\n");

	FreeArray2D(a);
	FreeArray2D(b);
	FreeArray2D(c);

	/*
	 * now test matrix multiply
	 */
	a = MakeArray2D(2,2);
	b = MakeArray2D(2,1);

	a->data[0*2+0] = 1.5;
	a->data[0*2+1] = 2.5;
	b->data[0] = 5.0;

	a->data[1*2+0] = 3.5;
	a->data[1*2+1] = 4.0;
	b->data[1] = 7.0;

	ac = CopyArray2D(a);
	c = RegressMatrix2D(a,b);
	d = MultiplyArray2D(ac,c);

	printf("a: ");
	PrintArray2D(a);
	printf("b: ");
	PrintArray2D(b);
	printf("c: ");
	PrintArray2D(c);
	printf("d: ");
	PrintArray2D(d);


	FreeArray2D(a);
	FreeArray2D(b);
	FreeArray2D(c);
	FreeArray2D(d);
	
	
	return(0);
}

#endif

	

	



