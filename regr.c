#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mioregress.h"

#define ARGS "x:y:"
char *Usage = "usage: regr -x xfile\n\
\t-y yfile\n";

char Xfile[4096];
char Yfile[4096];


int main(int argc, char *argv[])
{
	int i;
	int j;
	int c;
	MIO *d_mio;
	Array2D *x;
	MIO *xmio;
	Array2D *b;
	Array2D *y;
	MIO *ymio;
	unsigned long size;

	while((c = getopt(argc,argv,ARGS)) != EOF) {
		switch(c) {
			case 'x':
				strncpy(Xfile,optarg,sizeof(Xfile));
				break;
			case 'y':
				strncpy(Yfile,optarg,sizeof(Yfile));
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

	b = RegressMatrix2D(x,y);
	if(b == NULL) {
		fprintf(stderr,"regression failed\n");
		exit(1);
	}

	printf("b: ");
	PrintArray2D(b);

	printf("\n");

	FreeArray2D(x);
	FreeArray2D(y);
	FreeArray2D(b);

	return(0);
}


	

	



