CC=gcc
MPATH=../mio
CFLAGS=-g -I${MPATH}

LIBS=${MPATH}/mymalloc.o -lm

all: polyco-test polyco.o

polyco-test: polyco.c polyco.h
	${CC} ${CFLAGS} -DSTANDALONE -o polyco-test polyco.c ${LIBS}

polyco.o: polyco.c polyco.h
	${CC} ${CFLAGS} -c polyco.c

clean:
	rm *.o polyco-test 
