CC = g++
CFLAGS = -std=c++11 -Wall -pg -O0 -Wpedantic
SRC = msa.cpp matrix.cpp
HEADERS = msa.h matrix.h
TARGET = msa

$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f msa

.PHONY: clean


po