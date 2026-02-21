#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mioregress.h"

#define DEFAULT_RATE (0.01)
#define DEFAULT_MOMENTUM (0.01)
#define DEFAULT_ERR (0.1)

#define CLIP 1

/*
 * single-hidden layer version of multi-layer perceptron
 *
 * full interconnected
 *
 * hidden layer is same dimension as input
 *
 * single output layer (regression)
 *
 * https://www.cse.unsw.edu.au/~cs9417ml/MLP2/BackPropagation.html
 * http://www.dontveter.com/bpr/public2.html
 * https://medium.com/@tiago.tmleite/neural-networks-multilayer-perceptron-and-the-backpropagation-algorithm-a5cd5b904fde
 */

#define RAND() drand48()
#define SEED(x) srand48(x)

#define ARGS "x:y:TE:R:I:vt:"
char *Usage = "usage: mlp-regr -x xfile\n\
\t-y yfile\n\
\t-t testfile\n\
\t-T <training mode>\n\
\t-E error threshold (training mode)\n\
\t-R learning rate (training mode)\n\
\t-I max iterations\n\
\t-v verbose mode\n";

char Xfile[4096];
char Yfile[4096];
char Tfile[4096];
int Training;
double Error;
double Rate;
int Verbose;
int Iterations;

struct net_stc
{
	Array2D *x;
	Array2D *y;
	Array2D *ItoH; // w from input to hidden
	Array2D *biastoH; // bias hidden
	Array2D *HtoO; // w from hidden to output
	Array2D *biastoO; // bias output
	Array2D *Hx; // hidden outputs (including transfer)
	Array2D *Ox; // estimates
	double rate;
	double xmean;
	double xsd;
	double ymean;
	double ysd;
	double xmax;
	double xmin;
	double ymax;
	double ymin;
};
typedef struct net_stc Net;

void PrintNet(Net *n)
{
	Array2D *temp;
	printf("input vectors\n");
	PrintArray2D(n->x);
	printf("layer 1 weights\n");
	PrintArray2D(n->ItoH);
	printf("layer 1 biases\n");
	PrintArray2D(n->biastoH);
	temp = TransposeArray2D(n->Hx);
	printf("layer 1 coefficients\n");
	PrintArray2D(n->Hx);
	FreeArray2D(temp);
	printf("layer 2 weights\n");
	PrintArray2D(n->HtoO);
	printf("layer 2 biases\n");
	PrintArray2D(n->biastoO);
	temp = TransposeArray2D(n->Ox);
	printf("layer 2 coefficients\n");
	PrintArray2D(temp);
	FreeArray2D(temp);
	printf("layer 2 outputs\n");
	PrintArray2D(n->Ox);
	return;
}
	

	


void Randomize(Array2D *a)
{
	int i;
	int j;

	for(i=0; i < a->ydim; i++) {
		for(j=0; j < a->xdim; j++) {
		a->data[i*a->xdim+j] = (2*RAND()) -1 ;
//			a->data[i*a->xdim+j] = RAND();
			a->data[i*a->xdim+j] = (2*RAND()) - 1;
		}
	}

	return;
}

void ArrayMeanSD(Array2D *a, double *mu, double *sd)
{
	double mean;
	double var;
	double total;
	double acc;
	int i;
	int j;

	total = 0;
	acc = 0;
	for(i=0; i < a->ydim; i++) {
		for(j=0; j < a->xdim; j++) {
			acc += a->data[i*a->xdim+j];
			total++;
		}
	}
	mean = acc / total;

	total = 0;
	acc = 0;
	for(i=0; i < a->ydim; i++) {
		for(j=0; j < a->xdim; j++) {
			acc += ((a->data[i*a->xdim+j] - mean) * (a->data[i*a->xdim+j] - mean));
			total++;
		}
	}
	var = acc / total;

	*mu = mean;
	*sd = sqrt(var);
	return;
}

void ArrayMinMax(Array2D *a, double *min, double *max)
{
	double lmin = 99999999999999;
	double lmax = -99999999999999;
	int i;
	int j;

	for(i=0; i < a->ydim; i++) {
		for(j=0; j < a->xdim; j++) {
			if(a->data[i*a->xdim+j] < lmin) {
				lmin = a->data[i*a->xdim+j];
			}
			if(a->data[i*a->xdim+j] > lmax) {
				lmax = a->data[i*a->xdim+j];
			}
		}
	}

	*min = lmin;
	*max = lmax;
	return;
}

void ScaleArray2D(Array2D *a, double mean, double sd)
{
	int i;
	int j;
	for(i=0; i < a->ydim; i++) {
		for(j=0; j < a->xdim; j++) {
			a->data[i*a->xdim+j] = (a->data[i*a->xdim+j] - mean) / sd;
		}
	}
	return;
}

void ScaleArray2DMinMax(Array2D *a, double min, double max)
{
	int i;
	int j;
	for(i=0; i < a->ydim; i++) {
		for(j=0; j < a->xdim; j++) {
			a->data[i*a->xdim+j] = (a->data[i*a->xdim+j] - min) / (max - min);
		}
	}
	return;
}

void UnScaleArray2D(Array2D *a, double mean, double sd)
{
	int i;
	int j;
	for(i=0; i < a->ydim; i++) {
		for(j=0; j < a->xdim; j++) {
			a->data[i*a->xdim+j] = (a->data[i*a->xdim+j]*sd) + mean;
		}
	}
	return;
}

void UnScaleArray2DMinMax(Array2D *a, double min, double max)
{
	int i;
	int j;
	for(i=0; i < a->ydim; i++) {
		for(j=0; j < a->xdim; j++) {
			a->data[i*a->xdim+j] = (a->data[i*a->xdim+j]*(max-min)) + min;
		}
	}
	return;
}
Net *InitNet(Array2D *x, Array2D *y, double rate)
{
	Net *n;

	n = (Net *)malloc(sizeof(Net));
	if(n == NULL) {
		exit(1);
	}

	n->x = x;
	n->y = y;

	ArrayMeanSD(n->x,&n->xmean,&n->xsd);
	ScaleArray2D(n->x,n->xmean,n->xsd);
	ArrayMeanSD(n->y,&n->ymean,&n->ysd);
	ScaleArray2D(n->y,n->ymean,n->ysd);

#if 0
	ArrayMinMax(n->x,&n->xmean,&n->xsd);
	ScaleArray2DMinMax(n->x,n->xmean,n->xsd);
	ArrayMinMax(n->y,&n->ymean,&n->ysd);
	ScaleArray2DMinMax(n->y,n->ymean,n->ysd);
#endif

	n->ItoH = MakeArray2D(x->xdim,x->xdim);
	if(n->ItoH == NULL) {
		goto out;
	}
	Randomize(n->ItoH);

	n->HtoO = MakeArray2D(x->xdim,y->xdim);
	if(n->HtoO == NULL) {
		goto out;
	}
	Randomize(n->HtoO);

	n->biastoH = MakeArray1D(x->xdim);
	if(n->biastoH == NULL) {
		goto out;
	}
	Randomize(n->biastoH);

	n->biastoO = MakeArray1D(y->xdim);
	if(n->biastoO == NULL) {
		goto out;
	}
	Randomize(n->biastoO);

	n->Hx = MakeArray1D(x->xdim);
	if(n->Hx == NULL) {
		goto out;
	}
	n->Hx = MakeArray1D(x->xdim);
	if(n->Hx == NULL) {
		goto out;
	}
	n->Ox = MakeArray1D(y->xdim);
	if(n->Ox == NULL) {
		goto out;
	}

	n->rate = rate;

	return(n);
out:
	fprintf(stderr,"init could not create arrays\n");
	exit(1);
}

double GlobalErrorOld(Net *n, Array2D *yprime)
{
	int i;
	int j;
	double v1;
	double v2;

	double sum = 0;
	for(i=0; i < n->x->xdim; i++) {
		for(j=0; j < n->y->ydim; j++) {
			/*
			v1 = (n->y->data[j*n->y->xdim+i]*n->xsd)+n->xmean;
			v2 = (yprime->data[j*yprime->xdim + i]*n->xsd)+n->xmean;
			*/
			v1 = n->y->data[j*n->y->xdim+0];
			v2 = yprime->data[j*yprime->xdim + i];
printf("y: %f, yprime[%d %d]: %f\n",v1,j,i,v2);
			sum += ((v1 - v2) * (v1 - v2));
		}
	}
	return(sum/2.0);
}

double GlobalError(Net *n, Array2D *yprime)
{
	int i;
	int j;
	double v1;
	double v2;

	double sum = 0;
	for(j=0; j < n->y->ydim; j++) {
		/*
		v1 = (n->y->data[j*n->y->xdim+i]*n->xsd)+n->xmean;
		v2 = (yprime->data[j*yprime->xdim + i]*n->xsd)+n->xmean;
		*/
		v1 = n->y->data[j*n->y->xdim+0];
		v2 = yprime->data[j*yprime->xdim + 0];
//printf("y: %f, yprime[%d]: %f\n",v1,j,v2);
		sum += ((v1 - v2) * (v1 - v2));
	}
	return(sum/2.0);
}
void FeedForward(int input, 
		 Net *n,
		 Array2D *yprime)
				
{
	int i;
	int node;
	double sum;

	/*
 	 * for the hidden layer
 	 * for each node, compute the sum
 	 *
 	 * each row of x is a different input
 	 * each column of weights (ItoH and HtoO) are the weights for a node in a layer
 	 *
 	 */
	for(node = 0; node < n->ItoH->xdim; node++) {
		sum = 0;
		for(i=0; i < n->x->xdim; i++) {
			sum += (n->x->data[input*n->x->xdim + i] * n->ItoH->data[i*n->ItoH->xdim + node]);
		}
		/*
 		 * for sigmoid
		n->Hx->data[node] = Sigmoid(sum);
		 */
		n->Hx->data[node] = sum + n->biastoH->data[node];
	}


	/*
 	 * now the output layer
 	 * ItoH->xdim is the number of nodes in the hidden layer
 	 */
//	for(node = 0; node < n->HtoO->xdim; node++) {
	for(node = 0; node < n->ItoH->xdim; node++) {
		sum = 0;
		for(i=0; i < n->ItoH->xdim; i++) {
			sum += (n->Hx->data[i] * n->HtoO->data[i*n->HtoO->xdim + node]);
		}
		n->Ox->data[node] = sum + n->biastoO->data[node];
//printf("sum: %f, bias: %f\n",sum,n->biastoO->data[node]);
//printf("yprime: %d %d\n",input,node);
		yprime->data[input*yprime->xdim + node] = n->Ox->data[node];
	}

	return;
}

double BackPropagation(int input,
			Net *n)
{
	int i;
	int j;
	int node;
	int inode;
	int onode;
	double d;
	double delta;
	Array2D *out_delta;
	double sum;
	double de_dw;



	out_delta = MakeArray1D(n->HtoO->xdim);
	if(out_delta == NULL) {
		fprintf(stderr,"Backpropagation no space for carried deltas\n");
		exit(1);
	}
	/*
 	 * do the HtoO weights and bias first
 	 * Ox is the output
 	 * Hx is the computed output friom each hidden node
 	 */
	for(node=0; node < n->HtoO->xdim; node++) {
		delta = -2 * (n->y->data[node*n->y->xdim+input] - n->Ox->data[node]);
		/*
 		 * for sigmoid
		delta = d * n->Ox->data[node] * (1.0 - n->Ox->data[node]);
		 */
		//delta = d * n->Ox->data[node];
		if(delta > CLIP) {
			delta = CLIP;
		} else if (delta < -CLIP) {
			delta = -CLIP;
		}
		out_delta->data[node] = delta;
		for(inode=0; inode < n->ItoH->xdim; inode++) {
			de_dw = delta * n->Hx->data[inode];
			n->HtoO->data[inode*n->HtoO->xdim + node] -= (n->rate * de_dw);
		}
		n->biastoO->data[node] -= (n->rate * delta);
	}


	/*
 	 * now do the hidden layer
 	 */
	for(node=0; node < n->ItoH->xdim; node++) {
		sum = 0;
		for(onode=0; onode < n->HtoO->xdim; onode++) {
			sum += (n->HtoO->data[node*n->HtoO->xdim + onode] * out_delta->data[onode]);
		}
		delta = sum;
		/*
		 * for sigmoid
		d = (n->y->data[onode] - n->Ox->data[onode]) * n->Ox->data[onode] * sum;
		*/
		delta = sum;
		if(delta > CLIP) {
			delta = CLIP;
		} else if (delta < -CLIP) {
			delta = -CLIP;
		}
		for(inode=0; inode < n->x->xdim; inode++) {
			de_dw = delta * n->x->data[input*n->x->xdim + inode];
			n->ItoH->data[inode*n->ItoH->xdim + node] -= (n->rate * de_dw);
		}
		n->biastoH->data[node] -= n->rate * delta;
	}
	FreeArray2D(out_delta);

	return(0.0);
}

int main(int argc, char *argv[])
{
	int i;
	int j;
	int k;
	int c;
	int input;
	MIO *d_mio;
	Array2D *x;
	MIO *xmio;
	Array2D *y;
	MIO *ymio;
	Array2D *t;
	MIO *tmio;
	Array2D *yprime;
	double err;
	unsigned long size;
	Net *n;
	Array2D *temp;
	Array2D *temp1;
	int curr_iter;

	Iterations = 50000;

	while((c = getopt(argc,argv,ARGS)) != EOF) {
		switch(c) {
			case 'x':
				strncpy(Xfile,optarg,sizeof(Xfile));
				break;
			case 'y':
				strncpy(Yfile,optarg,sizeof(Yfile));
				break;
			case 'E':
				Error = atof(optarg);
				break;
			case 'R':
				Rate = atof(optarg);
				break;
			case 'T':
				Training = 1;
				break;
			case 'v':
				Verbose = 1;
				break;
			case 't':
				strncpy(Tfile,optarg,sizeof(Tfile));
				break;
			case 'I':
				Iterations = atoi(optarg);
				break;
			default:
				fprintf(stderr,
			"unrecognized command: %c\n",(char)c);
				fprintf(stderr,"%s",Usage);
				exit(1);
		}
	}

	if(Training == 1) {
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
	}

	if(Training == 1) {
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
	}

	if(Error == 0.0) {
		Error = DEFAULT_ERR;
	}
	if(Rate == 0.0) {
		Rate = DEFAULT_RATE;
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
	yprime = MakeArray2D(y->ydim,x->xdim);

	if(x == NULL) {
		fprintf(stderr,"no space for x array\n");
		exit(1);
	}

	if(y == NULL) {
		fprintf(stderr,"no space for y array\n");
		exit(1);
	}
	if(yprime == NULL) {
		fprintf(stderr,"no space for yprime array\n");
		exit(1);
	}

	if(Training == 1) {
		n = InitNet(x,y,Rate);

		err = 1000000;
		i = 0;
		curr_iter = 0;
		while((err > Error) && (curr_iter < Iterations)) {
			for(input=0; input < x->ydim; input++) {
				FeedForward(input,n,yprime);
//				if(Verbose) {
//					PrintNet(n);
//				}
				BackPropagation(input,n);
			}
			err = GlobalError(n,yprime);
			printf("iter: %d, err: %f %f\n",i,err,Error);
			if(Verbose) {
#if 0
				temp1=CopyArray2D(yprime);
				UnScaleArray2D(temp1,n->ymean,n->ysd);
				temp=CopyArray2D(y);
				UnScaleArray2D(temp,n->ymean,n->ysd);
				for(j=0; j < y->ydim; j++) {
					for(k=0; k < y->xdim; k++) {
						printf("%f %f (%f)\n",
							temp->data[j*y->xdim+k],
							temp1->data[j*yprime->xdim+k],
							temp->data[j*y->xdim+k] - temp1->data[j*yprime->xdim+k]);
					}
				}
				printf("\n");
				fflush(stdout);
				FreeArray2D(temp);
				FreeArray2D(temp1);
#endif
			}
			i++;
			curr_iter++;
		}
	}

	if(Tfile[0] != 0) {
		size = MIOSize(Tfile);
		d_mio = MIOOpenText(Tfile,"r",size);
		if(d_mio == NULL) {
			fprintf(stderr,"couldn't open %s\n",Tfile);
			exit(1);
		}
		tmio = MIODoubleFromText(d_mio,NULL);
		if(tmio == NULL) {
			fprintf(stderr,"no valid data in %s\n",Tfile);
			exit(1);
		}
		t = MakeArray2DFromMIO(tmio);
		if(t == NULL) {
			fprintf(stderr,"no space for test matrix\n");
			exit(1);
		}
		FreeArray2D(yprime);
		yprime = MakeArray2D(t->ydim,y->xdim);
		if(yprime == NULL) {
			fprintf(stderr,"no space for inference predictions\n");
			exit(1);
		}
		/*
 		 * make the test inputs the x inputs
 		 */
		FreeArray2D(n->x);
		x = t;
		n->x = t;
		ScaleArray2D(n->x,n->xmean,n->xsd);
		for(input=0; input < t->ydim; input++) {
			FeedForward(input,n,yprime);
		}
		UnScaleArray2D(yprime,n->ymean,n->ysd);
		for(i=0;i < yprime->ydim; i++) {
			for(j=0; j < n->y->xdim; j++) {
				printf("%f ",yprime->data[yprime->xdim*i+j]);
			}
			printf("\n");
		}
	}


	FreeArray2D(x);
	FreeArray2D(y);

	return(0);
}


