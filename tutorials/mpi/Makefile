CXX = g++
CXXFLAGS = -g -O3 -Wall -Wextra -Wpedantic -std=c++11

# needs its own set of compilers 
MPICXX?=/usr/local/bin/mpic++
MPILIB=-L/usr/local/lib/ -lmpi
MPIINCLUDE=-I/usr/local/include/

TARGETS= serial mpi

default : all

serial: serial.cpp
	$(CXX) $(CXXFLAGS) serial.cpp -o $@

mpi: mpi.cpp
	$(MPICXX) $(CXXFLAGS) $(MPILIB) $(MPIINCLUDE) mpi.cpp -o $@


all : $(TARGETS)

clean:
	rm -f $(TARGETS) *.o *.exe
