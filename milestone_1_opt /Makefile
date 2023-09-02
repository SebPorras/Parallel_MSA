CC = g++
CFLAGS = -std=c++11 -Wall -O3 -pg -Wpedantic -ffast-math
SRC = msaOpt.cpp matrixOpt.cpp
HEADERS = msaOpt.h matrixOpt.h
TARGET = msa_3

$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f msa_3
	
.PHONY: clean


