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
const int GAP = -10;

struct Sequences {
    int numSeqs;
    std::vector<std::string> seqs;
    std::map<std::string, int> ids; 
    std::map<int, std::string> idToName; 
};

void read_fasta_file(std::string fileName, std::unique_ptr<Sequences>& seqs); 
void calc_distances(std::vector<double>& distances, int matDims, 
        std::unique_ptr<Sequences>& seqs);
double perform_alignment(std::string seq1, std::string seq2);
void print_matrix(int* M, int rows, int cols);
double calculate_similarity(std::string seq1, std::string seq2); 
void print_lower_diagnonal(std::vector<double>& lowerD, int matDims); 

#endif
