
CC=g++
CFLAGS=-pedantic -Wall -Werror -std=gnu99 -g -pg

.PHONY: clean

all: msa

msa: msa.o 
	${CC} ${CLFAGS} msa.o -o msa

msa.o: msa.cpp msa.h
	${CC} ${CLFAGS} -o msa.o -c msa.cpp

clean: 
	rm -f *.o msa 
