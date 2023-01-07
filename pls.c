#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mioregress.h"
#include "redblack.h"

Array2D *PLS(Array2D *x, Array2D *y, int components)
{
	Array1D *p_k = NULL;
	Array1D *p0 = NULL;
	Array1D *t_k = NULL;
	Array1D *t_k1 = NULL;
	double tk;
	Array2D *x_k = NULL;
	Array2D *x_k1 = NULL;
	Array1D *w_k = NULL;
	Array1D *w_k1 = NULL;
	Array2D *P = NULL;
	Array2D *W = NULL;
	Array1D *q = NULL;
	double qk;
	Array1D *temp1 = NULL;
	Array2D *temp2 = NULL;
	Array1D *temp3 = NULL;
	Array2D *temp4 = NULL;
	Array1D *B = NULL;
	int k;
	int i;
	int j;
	double tempd;
	int tot_c;

	/*
	 * X_0 = x
	 */
	x_k = CopyArray2D(x);
	if(x_k == NULL) {
		goto out;
	}


	/*
	 * w_0 = x^t * y / ||x^t y||
	 */
	temp2 = TransposeArray2D(x_k);
	if(temp2 == NULL) {
		goto out;
	}
	/*
	 * temp1 is x^t * y
	 */
	temp1 = MultiplyArray2D(temp2,y);
	if(temp1 == NULL) {
		goto out;
	}
	w_k = MakeArray1D(temp1->ydim);
	if(w_k == NULL) {
		goto out;
	}

	tempd = NormArray1D(temp1);
	for(i=0; i < w_k->ydim; i++) {
		w_k->data[i*w_k->xdim+0] =
			temp1->data[i*temp1->xdim+0] / tempd;
	}

	t_k = MultiplyArray2D(x,w_k);
	if(t_k == NULL) {
		goto out;
	}

	FreeArray2D(temp1);
	temp1 = NULL;
	FreeArray2D(temp2);
	temp2 = NULL;

	W = MakeArray2D(x->xdim,components);
	if(W == NULL) {
		goto out;
	}

	P = MakeArray2D(x->xdim,components);
	if(P == NULL) {
		goto out;
	}

	q = MakeArray1D(components);
	if(q == NULL) {
		goto out;
	}

	tot_c = components;
	for(k=0; k < components; k++) {
		/*
		 * tk = t_k^t * t_k
		 */
		temp1 = TransposeArray2D(t_k);
		if(temp1 == NULL) {
			goto out;
		}
		temp3 = MultiplyArray2D(temp1,t_k);
		if(temp3 == NULL) {
			goto out;
		}
		tk = temp3->data[0];
		FreeArray1D(temp3);
		temp3 = NULL;
		FreeArray1D(temp1);
		temp1 = NULL;

		/*
		 * t_k = t_k / tk
		 */
		for(i=0; i < t_k->ydim; i++) {
			t_k->data[i*t_k->xdim+0] =
				t_k->data[i*t_k->xdim+0] / tk;
		}

		/*
		 * p_k = x_k^t * t_k
		 */
		temp2 = TransposeArray2D(x_k);
		if(temp2 == NULL) {
			goto out;
		}
		p_k = MultiplyArray2D(temp2,t_k);
		if(p_k == NULL) {
			goto out;
		}
		FreeArray2D(temp2);
		temp2 = NULL;

		/*
		 * save off p0 for the end
		 */
		if(k == 0) {
			p0 = CopyArray1D(p_k);
			if(p0 == NULL) {
				goto out;
			}
		}

		/*
		 * qk = y^t * t_k
		 */
		temp1 = TransposeArray2D(y);
		if(temp1 == NULL) {
			goto out;
		}
		temp3 = MultiplyArray2D(temp1,t_k);
		if(temp3 == NULL) {
			goto out;
		}
		qk = temp3->data[0];
		FreeArray2D(temp1);
		temp1 = NULL;
		FreeArray2D(temp3);
		temp3 = NULL;
		if(qk == 0) {
			tot_c = k;
			break;
		}

		/*
		 * x_k1 = x_k - tk*t_k*p_k^t
		 */
		temp1 = TransposeArray2D(p_k);
		if(temp1 == NULL) {
			goto out;
		}
		for(i=0; i < t_k->ydim; i++) {
			t_k->data[i*t_k->xdim+0] =
				t_k->data[i*t_k->xdim+0] * tk;
		}
		temp3 = MultiplyArray2D(t_k,temp1);
		if(temp3 == NULL) {
			goto out;
		}
		FreeArray1D(temp1);
		temp1 = NULL;
		x_k1 = SubtractArray2D(x_k,temp3);
		if(x_k1 == NULL) {
			goto out;
		}
		FreeArray2D(temp3);
		temp3 = NULL;

		/*
		 * w_k1 = x_k1^t * y
		 */
		temp2 = TransposeArray2D(x_k1);
		if(temp2 == NULL) {
			goto out;
		}
		w_k1 = MultiplyArray2D(temp2,y);
		if(w_k1 == NULL) {
			goto out;
		}
		FreeArray2D(temp2);
		temp2 = NULL;

		/*
		 * t_k1 = x_k1 * w_k1
		 */
		t_k1 = MultiplyArray2D(x_k1,w_k1);
		if(t_k1 == NULL) {
			goto out;
		}

		/*
		 * load up the arrays
		 */
		for(i=0; i < W->ydim; i++) {
			W->data[i*W->xdim+k] =
				w_k->data[i*w_k->xdim+0];
		}
		for(i=0; i < P->ydim; i++) {
			P->data[i*P->xdim+k] =
				p_k->data[i*p_k->xdim+0];
		}
		q->data[k*q->xdim+0] = qk;

		/*
		 * swap them
		 */
		FreeArray2D(x_k);
		x_k = x_k1;
		x_k1 = NULL;
		FreeArray2D(w_k);
		w_k = w_k1;
		w_k1 = NULL;
		FreeArray2D(t_k);
		t_k = t_k1;
		t_k1 = NULL;
	}

	B = MakeArray1D(x->xdim+1);
	if(B == NULL) {
		goto out;
	}

	/*
	 * B = W * (P^t * W)^-1*q
	 * for i > 0 and 
	 * B_0 = q_0 - P(0)^t * B
	 */
	temp2 = TransposeArray2D(P);
	if(temp2 == NULL) {
		FreeArray1D(B);
		B = NULL;
		goto out;
	}
	/*
	 * temp4 is P^t * W
	 */
	temp4 = MultiplyArray2D(temp2,W);
	if(temp4 == NULL) {
		FreeArray1D(B);
		B = NULL;
		goto out;
	}
	FreeArray2D(temp2);
	temp2 = NULL;

	/*
	 * temp2 is (P^t * W)^-1
	 */
	temp2 = InvertArray2D(temp4);
	if(temp2 == NULL) {
		FreeArray1D(B);
		B = NULL;
		goto out;
	}
	FreeArray2D(temp4);
	temp4 = NULL;


	/*
	 * temp4 is W * (P^t * W)^-1
	 */
	temp4 = MultiplyArray2D(W,temp2);
	if(temp4 == NULL) {
		FreeArray1D(B);
		B = NULL;
		goto out;
	}
	FreeArray2D(temp2);
	temp2 = NULL;

	/*
	 * temp2 = B = W * (P^t * W)^-1 * q
	 */
	temp2 = MultiplyArray2D(temp4,q);
	if(temp2 == NULL) {
		FreeArray1D(B);
		B = NULL;
		goto out;
	}
	FreeArray2D(temp4);
	temp4 = NULL;

	for(i=1; i < B->ydim; i++) {
		B->data[i*B->xdim+0] = temp2->data[(i-1)*temp2->xdim+0];
	}

	temp1 = TransposeArray2D(p0);
	if(temp1 == NULL) {
		FreeArray1D(B);
		B = NULL;
		goto out;
	}
	/*
	 * temp2 is B without the intercept yterm B_0
	 *
	 * temp3 = P(0)^t * B
	 */
	temp3 = MultiplyArray2D(temp1,temp2);
	if(temp3 == NULL) {
		FreeArray1D(B);
		B = NULL;
		goto out;
	}
	FreeArray1D(temp1);
	temp1 = NULL;
	FreeArray1D(temp2);
	temp2 = NULL;
	/*
	 * B0 = q0 - P(0)^t*B
	 */
	B->data[0*B->xdim+0] = q->data[0*q->xdim+0] -
		temp3->data[0];
	FreeArray1D(temp3);
	temp3 = NULL;

out:
	if(p_k != NULL) {
		FreeArray1D(p_k);
	}
	if(p0 != NULL) {
		FreeArray1D(p0);
	}
	if(t_k != NULL) {
		FreeArray1D(t_k);
	}
	if(t_k1 != NULL) {
		FreeArray1D(t_k1);
	}
	if(x_k != NULL) {
		FreeArray2D(x_k);
	}
	if(x_k1 != NULL) {
		FreeArray2D(x_k1);
	}
	if(w_k != NULL) {
		FreeArray1D(w_k);
	}
	if(w_k1 != NULL) {
		FreeArray1D(w_k1);
	}
	if(P != NULL) {
		FreeArray2D(P);
	}
	if(W != NULL) {
		FreeArray2D(W);
	}
	if(q != NULL) {
		FreeArray1D(q);
	}
	if(temp1 != NULL) {
		FreeArray1D(temp1);
	}
	if(temp2 != NULL) {
		FreeArray2D(temp2);
	}
	if(temp3 != NULL) {
		FreeArray1D(temp3);
	}
	if(temp4 != NULL) {
		FreeArray1D(temp4);
	}

	return(B);
}


		





	
	

	
	


#ifdef STANDALONE

#define ARGS "x:y:c:S"
char *Usage = "usage: pls -x xfile\n\
\t-c count <number of components to use>\n\
\t-y yfile\n";

char Xfile[4096];
char Yfile[4096];
int Components;
int Summary;

int main(int argc, char *argv[])
{
	int c;
	int size;
	MIO *d_mio;
	Array2D *x;
	Array2D *rx;
	Array2D *y;
	MIO *xmio;
	MIO *ymio;
	Array2D *b;
	double rsq;
	double rmse;
	int i;
	int j;

	Summary = 0;
	while((c = getopt(argc,argv,ARGS)) != EOF) {
		switch(c) {
			case 'x':
				strncpy(Xfile,optarg,sizeof(Xfile));
				break;
			case 'y':
				strncpy(Yfile,optarg,sizeof(Yfile));
				break;
			case 'c':
				Components = atoi(optarg);
				break;
			case 'S':
				Summary = 1;
				break;
			default:
				fprintf(stderr,
			"unrecognized command: %c\n",(char)c);
				fprintf(stderr,"%s",Usage);
				exit(1);
		}
	}

	if(Xfile[0] == 0) {
		fprintf(stderr,"must specify xfile\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
	}

	if(Yfile[0] == 0) {
		fprintf(stderr,"must specify yfile\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
	}

	size = MIOSize(Xfile);
	d_mio = MIOOpenText(Xfile,"r",size);
	if(d_mio == NULL) {
		fprintf(stderr,"couldn't open %s\n",Xfile);
		exit(1);
	}
	xmio = MIODoubleFromText(d_mio,NULL);
	if(xmio == NULL) {
		fprintf(stderr,"no valid data in %s\n",Xfile);
		exit(1);
	}

	size = MIOSize(Yfile);
	d_mio = MIOOpenText(Yfile,"r",size);
	if(d_mio == NULL) {
		fprintf(stderr,"couldn't open %s\n",Yfile);
		exit(1);
	}
	ymio = MIODoubleFromText(d_mio,NULL);
	if(ymio == NULL) {
		fprintf(stderr,"no valid data in %s\n",Yfile);
		exit(1);
	}

	x = MakeArray2DFromMIO(xmio);
	y = MakeArray2DFromMIO(ymio);

	if(Components == 0) {
		Components = x->xdim;
	}

	b = PLS(x,y,Components);

	/*
	 * b has y-intercept term
	 */
	rx = MakeArray2D(x->ydim,x->xdim+1);
	if(rx == NULL) {
		exit(1);
	}
	for(i=0; i < rx->ydim; i++) {
		rx->data[i*rx->xdim+0] = 1;
		for(j=1; j < rx->xdim; j++) {
			rx->data[i*rx->xdim+j] =
				x->data[i*x->xdim+(j-1)];
		}
	}

	if(Summary == 0) {
		PrintArray1D(b);
	}


        rsq = RSquared(rx,b,y);
	rmse = RMSE(rx,b,y);
	printf("R^2: %f RMSE: %f\n",rsq,rmse);

	FreeArray2D(rx);
	FreeArray2D(x);
	FreeArray2D(y);
	FreeArray2D(b);

	return(0);
}

#endif
