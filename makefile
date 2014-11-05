CC=gcc
MPATH=../mio
CFLAGS=-g -I${MPATH}

LIBS=${MPATH}/mymalloc.o ${MPATH}/mio.o -lm

all: polyco-test polyco.o mioarray-test

polyco-test: polyco.c polyco.h
	${CC} ${CFLAGS} -DSTANDALONE -o polyco-test polyco.c ${LIBS}

polyco.o: polyco.c polyco.h
	${CC} ${CFLAGS} -c polyco.c

mioarray-test: mioarray.c mioarray.h
	${CC} ${CFLAGS} -DSTANDALONE -o mioarray-test mioarray.c ${LIBS}

clean:
	rm *.o polyco-test 
