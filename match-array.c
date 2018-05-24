#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#include "mioarray.h"

#define MAXINTERVAL 1800


Array2D *ComputeMatchArray(Array2D *pred_series, Array2D *meas_series)
{
	int i;
	int j;
	int err;
	int mcount;
	int k;
	double p_ts;
	double m_ts;
	double v_ts;
	double next_ts;
	double pv;
	double mv;
	Array2D *matched_array;

	matched_array = MakeArray2D(pred_series->ydim,2);
	if(matched_array == NULL) {
		return(NULL);
	}

	i = 0;
	j = 0;
	k = 0;
	m_ts = meas_series->data[i*2+0];
	next_ts = meas_series->data[(i+1)*2+0];
	p_ts = pred_series->data[j*2+0];
	while((i+1) < meas_series->ydim) {
#if 0
		/*
		 * deal with drop out in the measured series
		 */
		if(fabs(next_ts-m_ts) > 1200) {
			/* assumes that next_ts is bigger than p_ts and m_ts is smaller */
			while(fabs(next_ts-p_ts) > 3700) {
				j++;
				if(j >= pred_series->ydim) {
					break;
				}
				p_ts = pred_series->data[j*2+0];
			}
		}
#endif
		if(fabs(p_ts-next_ts) < fabs(p_ts-m_ts)) {
			matched_array->data[k*2+0] = pred_series->data[j*2+1];
			matched_array->data[k*2+1] = meas_series->data[(i+1)*2+1];
			v_ts = next_ts;
		} else {
			matched_array->data[k*2+0] = pred_series->data[j*2+1];
			matched_array->data[k*2+1] = meas_series->data[i*2+1];
			v_ts = m_ts;
		}


		
		if(fabs(v_ts - p_ts) > MAXINTERVAL) {
			if(v_ts > p_ts) { // move predictions forward
				while(v_ts > p_ts) {
					j++;
					if(j >= pred_series->ydim) {
						break;
					}
					p_ts = pred_series->data[j*2+0];
				}
				if(j >= pred_series->ydim) {
					break;
				}
			} else {
				while(v_ts < p_ts) { // move meas forward
					i++;
					if((i+1) >= meas_series->ydim) {
						break;
					}
					m_ts = next_ts;
					next_ts = meas_series->data[(i+1)*2+0];
					v_ts = next_ts;
				}
				if((i+1) >= meas_series->ydim) {
					break;
				}
			}
			continue; /* go back and try again */
		}
/*
printf("MATCHED(%d): p: %10.10f %f m: %10.10f %f\n",
j,p_ts,matched_array->data[j*2+0],
v_ts,matched_array->data[j*2+1]);
fflush(stdout);
*/

		k++;
		i++;
		j++;
		if(j >= pred_series->ydim) {
/*
printf("MATCHED SHORT: j: %d, i: %d, k: %d pydim: %lu mydim: %lu\n",
j,i,k,pred_series->ydim,meas_series->ydim);
fflush(stdout);
*/
			matched_array->ydim = k;
			return(matched_array);
		}
/*
if(pred_series->data[j*2+0] < p_ts) {
fprintf(stderr,"PANIC: pred series out of order p: %10.0f next: %10.0f\n",
p_ts,pred_series->data[j*2+0]);
fflush(stderr);
FreeArray2D(matched_array);
return(NULL);
}
*/

		p_ts = pred_series->data[j*2+0];
		m_ts = meas_series->data[i*2+0];
		next_ts = meas_series->data[(i+1)*2+0];
		
		while(next_ts < p_ts) {
			m_ts = next_ts;
			i++;
			if((i+1) >= meas_series->ydim) {
				break;
			}
			next_ts = meas_series->data[(i+1)*2+0];
		}
		if((i+1) >= meas_series->ydim) {
			break;
		}
	}

/*
printf("MATCHED END: j: %d, i: %d, k: %d, pydim: %lu mydim: %lu\n",
j,i,k,pred_series->ydim,meas_series->ydim);
fflush(stdout);
*/
	if(j < pred_series->ydim) {
		matched_array->data[k*2+0] = pred_series->data[j*2+1];
		matched_array->data[k*2+1] = meas_series->data[i*2+1];
	}
	matched_array->ydim = k;


	return(matched_array);
			
}

#define ARGS "x:y:"
char *Usage = "usage: match-array -x xfile\n\
\t-y yfile\n";

char Xfile[4096];
char Yfile[4096];
double Confidence;


int main(int argc, char *argv[])
{
	int i;
	int j;
	int c;
	MIO *d_mio;
	Array2D *x;
	Array2D *cx;
	MIO *xmio;
	Array2D *b;
	Array2D *y;
	MIO *ymio;
	Array2D *f;
	unsigned long size;
	Array2D *m;

	Confidence = 0;
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

	m = ComputeMatchArray(x,y);

	if(m == NULL) {
		fprintf(stderr,"couldn't compute match array\n");
		exit(1);
	}

	for(i=0; i < m->ydim; i++) {
		printf("%10.10f %f\n",m->data[i*2+0],m->data[i*2+1]);
		fflush(stdout);
	}

	FreeArray2D(x);
	FreeArray2D(y);
	FreeArray2D(m);

	return(0);
}


	

	



