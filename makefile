CC=gcc
MPATH=../mio
APATH=.
LPATH=./lapack-3.5.0/lapacke/include/
CFLAGS=-g -I${MPATH} -I${APATH} -I${LPATH} -I/usr/local/include

LIBS=${MPATH}/mymalloc.o ${MPATH}/mio.o -lm
ALIB=${APATH}/mioarray.o

LLIB=-L./lapack-3.5.0 -L/usr/local/Cellar/gcc/4.9.2_1/lib/gcc/4.9/ -llapacke -llapack -lrefblas -ltmglib -lgfortran

all: polyco-test polyco.o mioregress.o mioarray-test mioregress-test regr

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

mioregress.o: mioregress.c mioregress.h mioarray.h
	${CC} ${CFLAGS} -DUSELAPACK -c mioregress.c
#	${CC} ${CFLAGS} -c mioregress.c

clean:
	rm *.o polyco-test mioregress-test mioarray-test
