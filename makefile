CC=gcc
MPATH=../mio
APATH=.
LPATH=./lapack-3.5.0/lapacke/include/
CFLAGS=-g -I${MPATH} -I${APATH} -I${LPATH}

LIBS=${MPATH}/mymalloc.o ${MPATH}/mio.o -lm
ALIB=${APATH}/mioarray.o

LLIB=./lapack-3.5.0/liblapacke.a

all: polyco-test polyco.o mioregress.o mioarray-test mioregress-test

polyco-test: polyco.c polyco.h
	${CC} ${CFLAGS} -DSTANDALONE -o polyco-test polyco.c ${LIBS}

polyco.o: polyco.c polyco.h
	${CC} ${CFLAGS} -c polyco.c

mioarray.o: mioarray.h mioarray.c
	${CC} ${CFLAGS} -c mioarray.c

mioarray-test: mioarray.c mioarray.h
	${CC} ${CFLAGS} -DSTANDALONE -o mioarray-test mioarray.c ${LIBS}

mioregress-test: mioregress.c mioregress.h ${APATH}/mioarray.h ${APATH}/mioarray.o
	${CC} ${CFLAGS} -DSTANDALONE -o mioregress-test mioregress.c ${ALIB} ${LIBS} ${LLIB}

mioregress.o: mioregress.c mioregress.h mioarray.h
	${CC} ${CFLAGS} -DUSELAPACK -c mioregress.c

clean:
	rm *.o polyco-test 
