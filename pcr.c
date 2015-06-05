#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mioregress.h"

Array2D *VarCoVarArray2D(Array2D *x)
{
	/*
	 * assumes that the variables x_i are in rows
	 */
	Array2D *vcv;
	double acc;
	double count;
	double mu_i;
	double mu_j;
	int i;
	int j;
	int k;
	double *xdata;
	int c;

	vcv = MakeArray2D(x->xdim,x->xdim);
	if(vcv == NULL) {
		return(NULL);
	}

	xdata = x->data;
	c = x->xdim;

	for(i=0; i < x->xdim; i++) {
		for(j=0; j < x->xdim; j++) {
			/*
			 * get mu from each column i and j
			 */
			mu_i = 0;
			mu_j = 0;
			count = 0;
			for(k=0; k < x->ydim; k++) {
				mu_i += xdata[k*c+i];
				mu_j += xdata[k*c+j];
				count++;
			}
			mu_i = mu_i / count;
			mu_j = mu_j / count;
			acc = 0;
			count = 0;
			for(k=0; k < x->ydim; k++) {
				acc += ((xdata[k*c+i] - mu_i) *
				        (xdata[k*c+j] - mu_j));
				count++;
			}
			vcv->data[i*x->xdim+j] = acc / count;
		}
	}

	return(vcv);
}

Array2D *PCArray2D(Array2D *x)
{
	Array2D *vcv;
	Array2D *eigen;
	Array2D *u;

	vcv = VarCoVarArray2D(x);
	if(vcv == NULL) {
		return(NULL);
	}
	eigen = EigenVectorArray2D(vcv);
	if(eigen == NULL) {
		FreeArray2D(vcv);
		return(NULL);
	}

	u = NormalizeColsArray2D(eigen);
	if(u == NULL) {
		FreeArray2D(vcv);
		FreeArray2D(eigen);
		return(NULL);
	}

	FreeArray2D(vcv);
	FreeArray2D(eigen);
	return(u);
}

#ifdef STANDALONE

#define ARGS "x:y:EC:"
char *Usage = "usage: pca -x xfile\n\
\t-C count <number of components to use>\n\
\t-E <explain variation>\n";

char Xfile[4096];
char Yfile[4096];
int Explain;
int Components;

double UnscaleB0(double y_bar, Array2D *beta, Array2D *cs)
{
	double b0;
	int i;

	b0 = y_bar;
	for(i=0; i < beta->ydim; i++) {
		b0 -= ((beta->data[i*beta->xdim+0] * cs->data[0*cs->xdim+i])) /
				cs->data[1*cs->xdim+i];
	}

	return(b0);
}

int main(int argc, char *argv[])
{
	int c;
	int size;
	MIO *d_mio;
	Array2D *x;
	Array2D *y;
	Array2D *xt;
	Array2D *xtx;
	MIO *xmio;
	MIO *ymio;
	Array2D *w;
	Array2D *w_c;
	Array2D *cw_c;
	Array2D *ev;
	Array2D *u;
	Array2D *u_c;
	Array2D *vcv;
	Array2D *gamma;
	Array2D *b;
	Array2D *b_star;
	Array2D *rx;
	Array2D *cs;
	Array2D *sx;
	double y_bar;
	double acc;
	double count;
	double frac;
	int i;
	int j;
	double rsq;
	double rmse;
	double b0;

	Explain = 0;
	while((c = getopt(argc,argv,ARGS)) != EOF) {
		switch(c) {
			case 'x':
				strncpy(Xfile,optarg,sizeof(Xfile));
				break;
			case 'y':
				strncpy(Yfile,optarg,sizeof(Yfile));
				break;
			case 'C':
				Components = atoi(optarg);
				break;
			case 'E':
				Explain = 1;
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

	/*
	 * get center/scale mean and s
	 */
	cs = CenterScale(x);
	if(cs == NULL) {
		exit(1);
	}

	/*
	 * scale and center the data
	 */
	sx = MakeArray2D(x->ydim,x->xdim);
	if(sx == NULL) {
		exit(1);
	}
	/*
	 * scale and center the data
	 */
	for(i=0; i < sx->ydim; i++) {
		for(j=0; j < sx->xdim; j++) {
			sx->data[i*sx->xdim+j] =
			(x->data[i*x->xdim+j] - cs->data[0*cs->xdim+j]) /
				cs->data[1*cs->xdim+j];
		}
	}

	u = PCArray2D(sx);
	if(u == NULL) {
		fprintf(stderr,"couldn't get principal components\n");
		exit(1);
	}

	/*
	 * measurements are in rows of x
	 */
	w = MultiplyArray2D(sx,u);
	if(y == NULL) {
		fprintf(stderr,"couldn't compute PC transform\n");
		exit(1);
	}

	/*
	for(i=0; i < w->ydim; i++) {
		for(j=0; j < w->xdim; j++) {
			printf("%10.10f ",w->data[i*w->xdim+j]);
		}
		printf("\n");
	}
	*/

	if(Explain == 1) {
		vcv = VarCoVarArray2D(sx);
		if(vcv == NULL) {
			fprintf(stderr,"couldn't get var-co var matrix\n");
			exit(1);
		}
		ev = EigenValueArray2D(vcv);
		if(ev == NULL) {
			fprintf(stderr,"no eigen values for xtx\n");
			exit(1);
		}
		printf("correlation matrix\n");
		PrintArray2D(vcv);
		printf("eigen vectors (in columns)\n");
		PrintArray2D(u);
		printf("eigen values\n");
		PrintArray1D(ev);
		printf("varfrac: ");
		acc = 0;
		for(i=0; i < ev->ydim; i++) {
			acc += ev->data[i*ev->xdim+0];
		}
		for(i=0; i < ev->ydim; i++) {
			frac = ev->data[i*ev->xdim+0] / acc;
			printf("%f ",frac);
		}
		printf("\n");
		FreeArray2D(ev);
		FreeArray2D(vcv);
	}

	/*
	 * make data subset for regression
	 */
	w_c = MakeArray2D(sx->ydim,Components);
	if(w_c == NULL) {
		exit(1);
	}

	for(i=0; i < w_c->ydim; i++) {
		for(j=0; j < w_c->xdim; j++) {
			w_c->data[i*w_c->xdim+j] =
			  w->data[i*w->xdim+j];
		}
	}

	/*
	 * RegressMatrix2D is destructive
	 */
	cw_c = CopyArray2D(w_c);
	if(cw_c == NULL) {
		exit(1);
	}
	gamma = RegressMatrix2D(cw_c,y);
        if(gamma == NULL) {
                fprintf(stderr,"regression failed\n");
                exit(1);
        }
	FreeArray2D(cw_c);

	/*
	 * create array #Component of u vectors
	 */
	u_c = MakeArray2D(u->ydim,Components);
	if(u_c == NULL) {
		exit(1);
	}
	for(i=0; i < u_c->ydim; i++) {
		for(j=0; j < u_c->xdim; j++) {
			u_c->data[i*u_c->xdim+j] = u->data[i*u->xdim+j];
		}
	}


	/*
	 * create PCR estimator from PCA eigen vectors
	 * (without intercept)
	 */
	b_star = MultiplyArray2D(u_c,gamma);
	if(b_star == NULL) {
		fprintf(stderr,"couldn't form PCR estimator\n");
		exit(1);
	}

	/*
	 * get y_bar
	*/
	y_bar = 0;
	count = 0;
	for(i=0; i < y->ydim; i++) {
		y_bar += y->data[i*y->xdim+0];
		count++;
	}
	y_bar = y_bar / count;

	/*
	 * estimate of intercept for unscaled data
	 */
	b0 = UnscaleB0(y_bar,b_star,cs);

	b = MakeArray1D(b_star->ydim+1);
	if(b == NULL) {
		exit(1);
	}
	b->data[0*b->xdim+0] = b0;
	for(i=1; i < b->ydim; i++) {
		b->data[i*b->xdim+0] = 
			b_star->data[(i-1)*b_star->xdim+0]
			 / cs->data[1*cs->xdim+(i-1)];
	}

	printf("b: (PCR estimator)\n");
	PrintArray2D(b);

	/*
	 * now make regrssion array to compute fit
	 */
	rx = MakeArray2D(x->ydim,x->xdim+1);
	if(rx == NULL) {
		exit(1);
	}
	for(i=0; i < rx->ydim; i++) {
		rx->data[i*rx->xdim+0] = 1;
		for(j=1; j < rx->xdim; j++) {
			rx->data[i*rx->xdim+j] = x->data[i*x->xdim+(j-1)];
		}
	}

        rsq = RSquared(rx,b,y);
	rmse = RMSE(rx,b,y);
	printf("R^2: %f RMSE: %f\n",rsq,rmse);

	

	FreeArray2D(u);
	FreeArray2D(u_c);
	FreeArray2D(x);
	FreeArray2D(w);
	FreeArray2D(w_c);
	FreeArray2D(y);
	FreeArray2D(b);
	FreeArray2D(gamma);
	FreeArray2D(cs);
	FreeArray2D(sx);
	FreeArray2D(rx);
	FreeArray2D(b_star);

	return(0);
}

#endif
