CC=gcc
MPATH=../mio
APATH=.
DPATH=../distributions
EPATH=../euca-cutils
#LPATH=./lapack-3.5.0/lapacke/include/
LPATH=./lapack-3.11.0/LAPACKE/include
CBPATH=./lapack-3.11.0/CBLAS/include
BPATH=./lapack-3.11.0/BLAS/include
#LPATH=/usr/local/opt/lapack/include
#BPATH=/usr/local/opt/openblas/include
#CFLAGS=-g -I${MPATH} -I${APATH} -I${LPATH} -I${EPATH} -I${DPATH} -I/usr/local/include -I/usr/local/opt/openblas/include/

LIBS=${MPATH}/mymalloc.o ${MPATH}/mio.o ${EPATH}/libutils.a -lm ${DPATH}/normal.o
ALIB=${APATH}/mioarray.o

#centos 7LLIB=-L./lapack-3.8.0 -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/ -L/usr/lib64/atlas -llapacke -llapack -ltatlas -lgfortran
#LLIB=-L./lapack-3.8.0 -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/ -L/usr/lib64/atlas -llapacke -llapack -ltatlas -lgfortran
#LLIB=-llapacke -llapack -lgfortran -L/usr/lib64 -L/usr/lib64/atlas -lsatlas -lblas -lm
#LLIB=-L./lapack-3.8.0 -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/ -L/usr/lib64/atlas -L/usr/lib64 -lcblas -lblas -llapacke -llapack -latlas -lgfortran
#LLIB=-L./lapack-3.8.0 -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/ -llapacke -llapack -lblas -lgfortran
#LLIB=-L./lapack-3.8.0 -L/usr/local/Cellar/gcc/8.1.0/lib/gcc/8/ -llapacke -llapack -lblas -lgfortran
#LLIB=-L./lapack-3.5.0 -L/usr/local/Cellar/gcc/6.1.0/lib/gcc/6/ -llapacke -llapack -lblas -lgfortran
#LLIB=-L./lapack-3.5.0 -L/usr/local/Cellar/gcc/4.9.2_1/lib/gcc/4.9/ -llapacke -llapack -lrefblas -lblas -ltmglib -lgfortran
#LLIB=-L./lapack-3.5.0 -L/usr/local/Cellar/gcc/5.3.0/lib/gcc/5/ -llapacke -llapack -lrefblas -lblas -ltmglib -lgfortran

# for OSX, bew install lapack, bew install openblas
LLIB=-L/opt/homebrew/opt/lapack/lib -L /opt/homebrew/opt/openblas/lib -llapacke -llapack -lopenblas
CFLAGS=-DUSELAPACK -g -I${LPATH} -I${BPATH} -I${CBPATH} -I${MPATH} -I${APATH} -I${BPATH} -I${EPATH} -I${DPATH} -I/usr/local/include

all: polyco-test polyco.o mioregress.o mioarray-test mioregress-test regr mioeigen-test pca pcr match-array mlp-regr slp-cat mlp-cat

polyco-test: polyco.c polyco.h
	${CC} ${CFLAGS} -DSTANDALONE -o polyco-test polyco.c ${LIBS}

polyco.o: polyco.c polyco.h
	${CC} ${CFLAGS} -c polyco.c

mioarray.o: mioarray.h mioarray.c
#	${CC} ${CFLAGS} -c mioarray.c
	${CC} ${CFLAGS} -DUSELAPACK -c mioarray.c

mioarray-test: mioarray.c mioarray.h
#	${CC} ${CFLAGS} -DSTANDALONE -o mioarray-test mioarray.c ${LIBS}
	${CC} ${CFLAGS} -DSTANDALONE -DUSELAPACK -o mioarray-test mioarray.c ${LIBS} ${LLIB}

mioregress-test: mioregress.c mioregress.h ${APATH}/mioarray.h ${APATH}/mioarray.o
#	${CC} ${CFLAGS} -DSTANDALONE -o mioregress-test mioregress.c ${ALIB} ${LIBS} ${LLIB}
	${CC} ${CFLAGS} -DSTANDALONE -DUSELAPACK -o mioregress-test mioregress.c ${ALIB} ${LIBS} ${LLIB}

regr: regr.c mioregress.h ${APATH}/mioarray.h ${APATH}/mioarray.o mioregress.o
#	${CC} ${CFLAGS} -o regr regr.c mioregress.o ${ALIB} ${LIBS} ${LLIB}
	${CC} ${CFLAGS} -DUSELAPACK -o regr regr.c mioregress.o ${ALIB} ${LIBS} ${LLIB}

mlp-regr: mlp-regr.c ${APATH}/mioarray.h ${APATH}/mioarray.o
#	${CC} ${CFLAGS} -o regr regr.c mioregress.o ${ALIB} ${LIBS} ${LLIB}
	${CC} ${CFLAGS} -DUSELAPACK -o mlp-regr mlp-regr.c ${ALIB} ${LIBS} ${LLIB}

slp-cat: slp-cat.c ${APATH}/mioarray.h ${APATH}/mioarray.o
#	${CC} ${CFLAGS} -o regr regr.c mioregress.o ${ALIB} ${LIBS} ${LLIB}
	${CC} ${CFLAGS} -DUSELAPACK -o slp-cat slp-cat.c ${ALIB} ${LIBS} ${LLIB}

mlp-cat: mlp-cat.c ${APATH}/mioarray.h ${APATH}/mioarray.o
#	${CC} ${CFLAGS} -o regr regr.c mioregress.o ${ALIB} ${LIBS} ${LLIB}
	${CC} ${CFLAGS} -DUSELAPACK -o mlp-cat mlp-cat.c ${ALIB} ${LIBS} ${LLIB}

match-array: match-array.c ${APATH}/mioarray.h ${APATH}/mioarray.o
#	${CC} ${CFLAGS} -o regr regr.c mioregress.o ${ALIB} ${LIBS} ${LLIB}
	${CC} ${CFLAGS} -o match-array match-array.c ${ALIB} ${LIBS} ${LLIB}

pca: pca.c ${APATH}/mioarray.h ${APATH}/mioarray.o
#	${CC} ${CFLAGS} -DSTANDALONE -o pca pca.c ${ALIB} ${LIBS} ${LLIB}
	${CC} ${CFLAGS} -DSTANDALONE -DUSELAPACK -o pca pca.c ${ALIB} ${LIBS} ${LLIB}

pcr: pcr.c ${APATH}/mioarray.h ${APATH}/mioarray.o mioregress.o mioregress.h
#	${CC} ${CFLAGS} -DSTANDALONE -o pcr pcr.c mioregress.o ${ALIB} ${LIBS} ${LLIB}
	${CC} ${CFLAGS} -DSTANDALONE -DUSELAPACK -o pcr pcr.c mioregress.o ${ALIB} ${LIBS} ${LLIB}

mioeigen-test: ${APATH}/mioarray.h ${APATH}/mioarray.o mioeigen-test.c
#	${CC} ${CFLAGS} -o mioeigen-test mioeigen-test.c ${ALIB} ${LIBS} ${LLIB}
	${CC} ${CFLAGS} -DUSELAPACK -o mioeigen-test mioeigen-test.c ${ALIB} ${LIBS} ${LLIB}

mioregress.o: mioregress.c mioregress.h mioarray.h
	${CC} ${CFLAGS} -DUSELAPACK -c mioregress.c
#	${CC} ${CFLAGS} -c mioregress.c

clean:
	rm *.o polyco-test mioregress-test mioarray-test
