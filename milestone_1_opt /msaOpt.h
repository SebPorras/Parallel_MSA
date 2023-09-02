#ifndef MSAOPT_H
#define MSAOPT_H

#include <string>
#include <vector> 
#include <map>
#include <memory>
#include <list> 
#include <cmath> 
#include <limits.h>
#include <iomanip>
#include <ios>
#include <iostream>
#include <float.h>
#include <fstream> 
#include <chrono>
using namespace std;
extern int blosum[20][20];

const int MAX_SEQ_LEN = 200;
const int FILENAME = 1;

const int NUM_LETTERS = 20; 
const int ROW_LEN = 24;
const int MATRIX_SIZE = 600; 
const int ASCII_OFFSET = -65; 

const int CLI_ERROR = 1;
const int FILE_ERROR = 2;

const int MATCH = 3;
const int MISMATCH = 0;
const int GAP = -3;

struct Sequence {
    std::string seq; //actual sequence 
    std::string id; //sequence name 
    int index; //where the sequence is in the matrix  
};

float mean_difference(std::vector<Sequence>& c1, std::vector<Sequence>& c2, 
        const int numPoints, vector<float> distanceMatrix); 
std::vector<Sequence> read_fasta_file(std::string fileName); 
void UPGMA(std::vector<std::vector<Sequence>>& clusters, 
        vector<float>& distanceMatrix, vector<int>& subMatrix);

void print_seqs(vector<vector<Sequence>> clusters); 
#endif
