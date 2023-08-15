#ifndef MATRIX_H 
#define MATRIX_H

#include "msa.h"
#include <fstream>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstddef>
#include <cstring>
#include <iterator>


void calc_distances(std::vector<double>& distances, int matDims, 
        std::unique_ptr<Sequences>& seqs);
double perform_alignment(std::string seq1, std::string seq2);
void print_matrix(std::vector<double>& M, int dims);
double calculate_similarity(std::string seq1, std::string seq2); 
void print_lower_diagnonal(std::vector<double>& lowerD, int matDims); 
std::vector<int> create_matrix(std::string seq1, std::string seq2);
void align_seqs( std::string& seq1, std::string& seq2, std::string& aSeq1, 
        std::string& aSeq2, std::vector<int>& M, int rows, int cols);
int max_index(std::vector<double>& distances, int matDims); 

#endif 
