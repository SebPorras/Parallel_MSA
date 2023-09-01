CC = g++
CFLAGS = -std=c++11 -Wall -O0 -pg -Wpedantic -ffast-math
SRC = msaOpt.cpp matrixOpt.cpp
HEADERS = msaOpt.h matrixOpt.h
TARGET = msaOpt

$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f msaOpt
	
.PHONY: clean


