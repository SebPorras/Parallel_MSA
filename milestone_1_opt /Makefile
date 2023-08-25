CC = g++
CFLAGS = -std=c++11 -Wall -pg -O0 -Wpedantic
SRC = msaOpt.cpp matrixOpt.cpp
HEADERS = msaOpt.h matrixOpt.h
TARGET = msaOpt

$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f msaOpt
	
.PHONY: clean


