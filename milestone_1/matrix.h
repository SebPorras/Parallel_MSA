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
void print_matrix(int* M, int rows, int cols);
double calculate_similarity(std::string seq1, std::string seq2); 
void print_lower_diagnonal(std::vector<double>& lowerD, int matDims); 

#endif 
