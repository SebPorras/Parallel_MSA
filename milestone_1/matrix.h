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


double perform_alignment(std::string seq1, std::string seq2);
double calculate_similarity(std::string seq1, std::string seq2); 
std::vector<int> create_matrix(std::string seq1, std::string seq2);
void align_seqs( std::string& seq1, std::string& seq2, std::string& aSeq1, 
        std::string& aSeq2, std::vector<int>& M, int rows, int cols);
void calc_distances(int numSeqs, std::vector<Sequence>& seqs);
#endif 
