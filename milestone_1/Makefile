
CC=gcc
CFLAGS=-pedantic -Wall -Werror -std=gnu99 -g -pg

.PHONY: clean

all: msa

msa: msa.o io.o

io.o: io.h

clean: 
	rm -f *.o msa 