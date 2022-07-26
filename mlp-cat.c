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
 * https://medium.com/@tiago.tmleite/neural-networks-multilayer-perceptron-and-the-backpropagation-algorithm-a5cd5b904fde
 */

#define RAND() (drand48())
#define SEED(x) srand48(x)

#define ARGS "x:y:TE:R:M:vt:H:"
char *Usage = "usage: mlp-regr -x xfile\n\
\t-y yfile\n\
\t-t testfile\n\
\t-T <training mode>\n\
\t-H hidden_layers\n\
\t-E error threshold (training mode)\n\
\t-R learning rate (training mode)\n\
\t-M momentum (training mode)\n\
\t-v verbose mode\n";

char Xfile[4096];
char Yfile[4096];
char Tfile[4096];
int Training;
double Error;
double Rate;
double Momentum;
int Verbose;
int Hidden;

struct layer_stc
{
	Array2D *ItoL;		// array of weights on inputs to hidden layer
	Array2D *biastoL;	// biases to layer
	Array2D *Ox;		// layer outputs
};

typedef struct layer_stc LAYER;
	

struct net_stc
{
	Array2D *x;
	Array2D *y;
	Array1D *Ox;	// final output layer
	Array2D *HtoO;  // final weights
	Array1D *biastoO;  // final biases
	LAYER *layers;  // vector of LAYER structs, 1 for each hidden layer
	int layer_count;
	
	double rate;
	double mu;
	double xmean;
	double xsd;
	double ymean;
	double ysd;
};
typedef struct net_stc Net;

void PrintLayer(LAYER *layer)
{
	Array2D *temp;
	printf("weights\n");
	PrintArray2D(layer->ItoL);
	printf("biases\n");
	PrintArray2D(layer->biastoL);
	printf("outputs\n");
	temp = TransposeArray2D(layer->Ox);
	PrintArray2D(temp);
	FreeArray2D(temp);
	return;
}

void PrintNet(Net *n)
{
	Array2D *temp;
	int i;

	printf("input vectors\n");
	PrintArray2D(n->x);

	for(i=0; i < n->layer_count; i++) {
		printf("hidden layer %d\n",i);
		PrintLayer(&(n->layers[i]));
	}

	printf("final layer\n");
	printf("weights\n");
	PrintArray2D(n->HtoO);
	printf("biases\n");
	PrintArray2D(n->biastoO);
	printf("outputs\n");
	PrintArray2D(n->Ox);
	
	return;
}
	

double Sigmoid(double x)
{
	return(1.0 / (1.0 + exp(-1.0*x)));
}

void Randomize(Array2D *a)
{
	int i;
	int j;

	for(i=0; i < a->ydim; i++) {
		for(j=0; j < a->xdim; j++) {
			a->data[i*a->xdim+j] = (RAND()*2.0)-1.0;
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

Net *InitNet(Array2D *x, Array2D *y, int layers, double rate, double momentum)
{
	Net *n;
	int i;

	n = (Net *)malloc(sizeof(Net));
	if(n == NULL) {
		goto out;
	}

	n->layers = (LAYER *)malloc(layers * sizeof(LAYER));
	if(n->layers == NULL) {
		goto out;
	}
	n->layer_count = layers;
	for(i=0; i < n->layer_count; i++) {
		n->layers[i].ItoL = MakeArray2D(x->xdim,x->xdim);
		if(n->layers[i].ItoL == NULL) {
			goto out;
		}
		Randomize(n->layers[i].ItoL);
		n->layers[i].biastoL = MakeArray1D(x->xdim);
		if(n->layers[i].biastoL == NULL) {
			goto out;
		}
		Randomize(n->layers[i].biastoL);
		n->layers[i].Ox = MakeArray1D(x->xdim);
		if(n->layers[i].Ox == NULL) {
			goto out;
		}
	}

	n->x = x;
	n->y = y;

	ArrayMeanSD(n->x,&n->xmean,&n->xsd);
	ScaleArray2D(n->x,n->xmean,n->xsd);
	//ArrayMeanSD(n->y,&n->ymean,&n->ysd);
	//ScaleArray2D(n->y,n->ymean,n->ysd);
	//ScaleArray2D(n->y,n->xmean,n->xsd);

	n->Ox = MakeArray1D(y->xdim);
	if(n->Ox == NULL) {
		goto out;
	}

	n->HtoO = MakeArray2D(x->xdim,y->xdim);
        if(n->HtoO == NULL) {
                goto out;
        }
        Randomize(n->HtoO);

	n->biastoO = MakeArray1D(y->xdim);
        if(n->biastoO == NULL) {
                goto out;
        }
	Randomize(n->biastoO);

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
	double v1;
	double v2;

	double sum = 0;
	for(i=0; i < n->y->ydim; i++) {
		for(j=0; j < n->y->ydim; j++) {
			/*
			v1 = (n->y->data[j*n->y->xdim+i]*n->xsd)+n->xmean;
			v2 = (yprime->data[j*yprime->xdim + i]*n->xsd)+n->xmean;
			*/
			v1 = n->y->data[j*n->y->xdim+i];
			v2 = yprime->data[j*yprime->xdim + i];

			/*
 			 * assume p(1) is p(success) and p(-1) is p(failure) 
 			 */
			if(v1 == -1) {
				v1 = 0;
			}
			sum += ((v1 - v2) * (v1 - v2));
		}
	}
	return(sum/2.0);
}



void FeedForward(int input, 
		 Net *n,
		 Array2D *yprime)
				
{
	int i;
	int l;
	int node;
	double sum;

	/*
 	 * for the hidden layers
 	 * for each node, compute the sum
 	 *
 	 * each row of x is a different input
 	 * each column of weights (ItoH and HtoO) are the weights for a node in a layer
 	 *
 	 */
	for(l=0; l < n->layer_count; l++) {
		for(node = 0; node < n->layers[l].ItoL->xdim; node++) {
			sum = 0;
			if(l == 0) {
				for(i=0; i < n->x->xdim; i++) {
					sum += (n->x->data[input*n->x->xdim + i] * n->layers[l].ItoL->data[i*n->layers[l].ItoL->xdim + node]);
				} 
			} else {
				for(i=0; i < n->layers[l-1].Ox->xdim; i++) {
					sum += (n->layers[l-1].Ox->data[i] * n->layers[l].ItoL->data[i*n->layers[l].ItoL->xdim + node]);
				} 
			}
			/*
			 * for sigmoid
			 */
			n->layers[l].Ox->data[node] = Sigmoid(sum + n->layers[l].biastoL->data[node]);
		}
	}


	/*
 	 * now the output layer
 	 * ItoL->xdim is the number of nodes in the hidden layer
 	 */
	l = n->layer_count - 1;
	for(node = 0; node < n->HtoO->xdim; node++) {
		sum = 0;
		for(i=0; i < n->layers[l].ItoL->xdim; i++) {
			sum += (n->layers[l].Ox->data[i] * n->HtoO->data[i*n->HtoO->xdim + node]);
		}
		n->Ox->data[node] = Sigmoid(sum + n->biastoO->data[node]);
		yprime->data[input*yprime->xdim + node] = n->Ox->data[node];
	}

	return;
}

double BackPropagation(int input,
			Net *n)
{
	int i;
	int j;
	int l;
	int node;
	int inode;
	int onode;
	double d;
	double delta;
	Array2D *out_delta;
	Array2D *out_weights;
	Array2D *new_delta;
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
	l = n->layer_count - 1;
	for(node=0; node < n->HtoO->xdim; node++) {
		delta = -2 * (n->y->data[node*n->y->xdim+input] - n->Ox->data[node]) *
				n->Ox->data[node] * (1 - n->Ox->data[node]);
		/*
 		 * for sigmoid
		delta = d * n->Ox->data[node] * (1.0 - n->Ox->data[node]);
		 */
		if(delta > CLIP) {
			delta = CLIP;
		} else if (delta < -CLIP) {
			delta = -CLIP;
		}
		out_delta->data[node] = delta;
		for(inode=0; inode < n->layers[l].ItoL->xdim; inode++) {
			de_dw = delta * n->layers[l].Ox->data[inode];
			n->HtoO->data[inode*n->HtoO->xdim + node] -= (n->rate * de_dw);
		}
		n->biastoO->data[node] -= (n->rate * delta);
	}


	/*
 	 * now do the hidden layers
 	 */
	out_weights = n->HtoO;
	for(l=n->layer_count-1; l >=0; l--) {
		new_delta = MakeArray1D(n->layers[l].ItoL->xdim);
		if(new_delta == NULL) {
			exit(1);
		}
		for(node=0; node < n->layers[l].ItoL->xdim; node++) {
			sum = 0;
			for(onode=0; onode < out_weights->xdim; onode++) {
				sum += (out_weights->data[node*out_weights->xdim + onode] * out_delta->data[onode]);
			}
			delta = (sum * n->layers[l].Ox->data[node] * (1 - n->layers[l].Ox->data[node]));
			new_delta->data[node] = delta;
			 /*
			 * for sigmoid
			d = (n->y->data[onode] - n->Ox->data[onode]) * n->Ox->data[onode] * sum;
			*/
			if(delta > CLIP) {
				delta = CLIP;
			} else if (delta < -CLIP) {
				delta = -CLIP;
			}
			for(inode=0; inode < n->x->xdim; inode++) {
				if(l == 0) {
					de_dw = delta * n->x->data[input*n->x->xdim + inode];
				} else {
					de_dw = delta * n->layers[l-1].Ox->data[inode];
				}
				n->layers[l].ItoL->data[inode*n->layers[l].ItoL->xdim + node] -= (n->rate * de_dw);
			}
			n->layers[l].biastoL->data[node] -= n->rate * delta;
		}
		out_weights = n->layers[l].ItoL;
		FreeArray2D(out_delta);
		out_delta = new_delta;
	}
	FreeArray2D(out_delta);

	return;
}

int main(int argc, char *argv[])
{
	int i;
	int j;
	int k;
	int c;
	int iter;
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

	Hidden = 1;
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
			case 'H':
				Hidden = atoi(optarg);
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
			case 't':
				strncpy(Tfile,optarg,sizeof(Tfile));
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
		n = InitNet(x,y,Hidden,Rate,Momentum);

		err = 100000;
		i = 0;
		iter = 0;
		while((err > Error) && (iter < 1000000)) {
			for(input=0; input < x->ydim; input++) {
				FeedForward(input,n,yprime);
				/*
				if(Verbose) {
					PrintNet(n);
				}
				*/
				BackPropagation(input,n);
			}
			err = GlobalError(n,yprime);
			if(Verbose) {
				temp1=CopyArray2D(yprime);
				//UnScaleArray2D(temp1,n->xmean,n->xsd);
				temp=CopyArray2D(y);
				//UnScaleArray2D(temp,n->xmean,n->xsd);
				printf("iter: %d, err: %f %f\n",i,err,Error);
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
			}
			i++;
			iter++;
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
		//UnScaleArray2D(yprime,n->xmean,n->xsd);
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


