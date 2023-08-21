
CC=g++
CFLAGS=-pedantic -Wall -Werror -std=gnu99 -g -pg -gp 

.PHONY: clean

all: msa

msa: msa.o matrix.o
	${CC} ${CLFAGS} msa.o matrix.o -o msa

msa.o: msa.cpp msa.h
	${CC} ${CLFAGS} -o msa.o -c msa.cpp

matrix.o: matrix.cpp matrix.h
	${CC} ${CLFAGS} -o matrix.o -c matrix.cpp


clean: 
	rm -f *.o msa 
