#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mioregress.h"

#define DEFAULT_RATE (0.01)
#define DEFAULT_MOMENTUM (0.1)
#define DEFAULT_ERR (0.1)

#define CLIP 10

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
 */

#define RAND() drand48()
#define SEED(x) srand48(x)

#define ARGS "x:y:TE:R:M:v"
char *Usage = "usage: mlp-regr -x xfile\n\
\t-y yfile\n\
\t-T <training mode>\n\
\t-E error threshold (training mode)\n\
\t-R learning rate (training mode)\n\
\t-M momentum (training mode)\n\
\t-v verbose mode\n";

char Xfile[4096];
char Yfile[4096];
int Training;
double Error;
double Rate;
double Momentum;
int Verbose;

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
	double mu;
	Array2D *prevItoH; // prev delta for back prop
	Array2D *prevbiastoH; // prev delta for back prop
	Array2D *prevHtoO; // prev delta for backprop
	Array2D *prevbiastoO; // prev delta for backprop
};
typedef struct net_stc Net;

void Randomize(Array2D *a)
{
	int i;
	int j;

	for(i=0; i < a->ydim; i++) {
		for(j=0; j < a->xdim; j++) {
			a->data[i*a->xdim+j] = RAND();
		}
	}

	return;
}

Net *InitNet(Array2D *x, Array2D *y, double rate, double momentum)
{
	Net *n;

	n = (Net *)malloc(sizeof(Net));
	if(n == NULL) {
		exit(1);
	}

	n->x = x;
	n->y = y;

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

	n->prevItoH = MakeArray2D(x->xdim,x->xdim);
	if(n->prevItoH == NULL) {
		goto out;
	}
	Randomize(n->prevItoH);

	n->prevHtoO = MakeArray2D(x->xdim,y->xdim);
	if(n->prevHtoO == NULL) {
		goto out;
	}
	Randomize(n->prevHtoO);

	n->prevbiastoH = MakeArray1D(x->xdim);
	if(n->prevbiastoH == NULL) {
		goto out;
	}
	Randomize(n->prevbiastoH);

	n->prevbiastoO = MakeArray1D(y->xdim);
	if(n->prevbiastoO == NULL) {
		goto out;
	}
	Randomize(n->prevbiastoO);

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
	n->mu = momentum;

	return(n);
out:
	fprintf(stderr,"init could not create arrays\n");
	exit(1);
}

double GlobalError(Net *n, Array2D *yprime)
{
	int i;
	int j;
	double sum = 0;
	for(i=0; i < n->y->ydim; i++) {
		for(j=0; j < n->y->ydim; j++) {
			sum += ((n->y->data[j*n->y->xdim+i] - yprime->data[j*yprime->xdim + i]) *
				(n->y->data[j*n->y->xdim+i] - yprime->data[j*yprime->xdim + i]));
		}
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
	for(node = 0; node < n->HtoO->xdim; node++) {
		sum = 0;
		for(i=0; i < n->ItoH->xdim; i++) {
			sum += (n->Hx->data[i] * n->HtoO->data[i*n->HtoO->xdim + node]);
		}
		n->Ox->data[node] = sum + n->biastoO->data[node];
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
	Array2D *curr_delta;
	Array1D *curr_bias;
	double sum;



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
	curr_delta = MakeArray2D(n->HtoO->ydim,n->HtoO->xdim);
	if(curr_delta == NULL) {
		fprintf(stderr,"Backpropagation no space for output deltas\n");
		exit(1);
	}
	curr_bias = MakeArray1D(n->biastoO->ydim);
	if(curr_bias == NULL) {
		fprintf(stderr,"Backpropagation no space for output bias\n");
		exit(1);
	}
	for(node=0; node < n->HtoO->xdim; node++) {
		d = (n->y->data[node*n->y->xdim+input] - n->Ox->data[node]);
		/*
 		 * for sigmoid
		delta = d * n->Ox->data[node] * (1.0 - n->Ox->data[node]);
		 */
		//delta = d * n->Ox->data[node];
		delta = d;
		if(delta > CLIP) {
			delta = CLIP;
		} else if (delta < -CLIP) {
			delta = -CLIP;
		}
		out_delta->data[node] = delta;
		for(inode=0; inode < n->ItoH->xdim; inode++) {
			n->HtoO->data[inode*n->HtoO->xdim + node] +=
				((n->rate * delta * n->Hx->data[inode]) + (n->mu*n->prevHtoO->data[inode*n->HtoO->xdim + node]));
			curr_delta->data[inode*n->HtoO->xdim + node] = n->rate * delta * n->Hx->data[inode] + 
					n->mu*n->prevHtoO->data[inode*n->HtoO->xdim + node];
		}
		n->biastoO->data[node] += ((n->rate * delta) + (n->mu*n->prevbiastoO->data[node]));
		curr_bias->data[node] = (n->rate * delta) + (n->mu*n->prevbiastoO->data[node]);
	}

	FreeArray2D(n->prevHtoO);
	n->prevHtoO = curr_delta;
	FreeArray2D(n->prevbiastoO);
	n->prevbiastoO = curr_bias;

	/*
 	 * now do the hidden layer
 	 */
	curr_delta = MakeArray2D(n->ItoH->ydim,n->ItoH->xdim);
	if(curr_delta == NULL) {
		fprintf(stderr,"Backpropagation no space for hidden deltas\n");
		exit(1);
	}
	curr_bias = MakeArray1D(n->biastoO->ydim);
	if(curr_bias == NULL) {
		fprintf(stderr,"Backpropagation no space for hidden bias\n");
		exit(1);
	}
	for(node=0; node < n->ItoH->xdim; node++) {
		sum = 0;
		for(onode=0; onode < n->HtoO->xdim; onode++) {
			sum += (n->HtoO->data[node*n->HtoO->xdim + onode] * out_delta->data[onode]);
		}
		/*
		 * for sigmoid
		d = (n->y->data[onode] - n->Ox->data[onode]) * n->Ox->data[onode] * sum;
		*/
		for(onode=0; onode < n->HtoO->xdim; onode++) {
			// d = (n->y->data[input*n->y->xdim + onode] - n->Ox->data[onode]) * sum;
			// delta = n->rate * d * n->x->data[input*n->x->xdim + node];
			delta = sum;
			if(delta > CLIP) {
				delta = CLIP;
			} else if (delta < -CLIP) {
				delta = -CLIP;
			}
			for(inode=0; inode < n->x->xdim; inode++) {
				n->ItoH->data[inode*n->ItoH->xdim + node] += 
						((n->rate*delta*n->x->data[input*n->x->xdim + node]) + 
							(n->prevItoH->data[inode*n->ItoH->xdim + node] * n->mu));
				curr_delta->data[inode*n->ItoH->xdim + node] = delta; 
			}
			//d = (n->y->data[input*n->y->xdim + onode ] - n->Ox->data[onode]) * n->Ox->data[onode] * n->biastoH->data[node];
			// d = n->biastoH->data[node] * out_delta->data[onode];
			delta = n->biastoO->data[onode] * out_delta->data[onode];
			if(delta > CLIP) {
				delta = CLIP;
			} else if (delta < -CLIP) {
				delta = -CLIP;
			}
			n->biastoH->data[node] += ((n->rate*delta) + (n->prevbiastoH->data[node] * n->mu));
			curr_bias->data[node] = (n->rate*delta) + (n->prevbiastoH->data[node] * n->mu);
		}
	}
	FreeArray2D(out_delta);
	FreeArray2D(n->prevItoH);
	n->prevItoH = curr_delta;
	FreeArray2D(n->prevbiastoH);
	n->prevbiastoH = curr_bias;

	return;
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
	Array2D *yprime;
	double err;
	unsigned long size;
	Net *n;

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
			case 'M':
				Momentum = atof(optarg);
				break;
			case 'T':
				Training = 1;
				break;
			case 'v':
				Verbose = 1;
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
	}
	if(Yfile[0] == 0) {
		fprintf(stderr,"must specify yfile\n");
		fprintf(stderr,"%s",Usage);
		exit(1);
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
	if(Momentum == 0.0) {
		Momentum = DEFAULT_MOMENTUM;
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
	yprime = MakeArray2D(y->ydim,y->xdim);

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
		n = InitNet(x,y,Rate,Momentum);

		err = 1000000;
		i = 0;
		while(err > Error) {
			for(input=0; input < x->ydim; input++) {
				FeedForward(input,n,yprime);
				BackPropagation(input,n);
			}
			err = GlobalError(n,yprime);
			if(Verbose) {
				printf("iter: %d, err: %f\n",i,err);
				for(j=0; j < y->ydim; j++) {
					for(k=0; k < y->xdim; k++) {
						printf("%f %f (%f)\n",
							y->data[j*y->xdim+k],
							yprime->data[j*yprime->xdim+k],
							y->data[j*y->xdim+k] - yprime->data[j*yprime->xdim+k]);
					}
				}
				printf("\n");
				fflush(stdout);
			}
			i++;
		}
	}


	FreeArray2D(x);
	FreeArray2D(y);

	return(0);
}


