#ifndef MSA_H
#define MSA_H

#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <algorithm>
#include <map>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector> 

const int MAX_SEQ_LEN = 200;
const int FILENAME = 1;

const int CLI_ERROR = 1;
const int FILE_ERROR = 2;


const int MATCH = 3;
const int MISMATCH = 0;
const int PENALTY = -3;
const int GAP = -3;

struct Sequences {
    int numSeqs;
    std::vector<std::string> seqs;
    std::map<std::string, int> ids; 
};


void read_fasta_file(std::string fileName, std::unique_ptr<Sequences>& seqs); 
void calc_dist(int* distances, int i, int j, std::unique_ptr<Sequences>& seqs);
void perform_alignment(std::string seq1, std::string seq2);
void print_matrix(int* M, int rows, int cols);
float calculate_similarity(std::string seq1, std::string seq2); 

#endif
