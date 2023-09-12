CXX?=g++
CXXFLAGS?=-g -O3 -Wall -Wextra -Wpedantic -std=c++11

NVCC?=nvcc
NVFLAGS?=-g -arch=sm_80 -Wno-deprecated-gpu-targets -O3

TARGETS=serial cuda_v1 cuda_v2 cuda_v3 cuda_v4

default : all

%.o : %.cpp
	$(XX) $(CXXFLAGS) $(CXXFLAGS_$(basename $<)) -c $< -o $@

% : %.cpp
	$(CXX) $(CXXFLAGS) $(CXXFLAGS_$@) $(filter %.o %.cpp, $^) $(LDFLAGS) $(LIBS_$@) $(LIB) -o $@

%.o : %.cu
	$(NVCC) $(NVFLAGS) $(NVFLAGS_$(basename $<)) -c $< -o $@

% : %.cu
	$(NVCC) $(NVFLAGS) $(NVFLAGS_$@) $(filter %.o %.cu, $^) $(LDFLAGS) $(LIBS_$@) $(LIB) -o $@

all : $(TARGETS)

clean:
	rm -f $(TARGETS) *.o
