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
#include <unordered_map>


double run_pairwise_alignment(Sequence& seq1, Sequence& seq2, bool modify);

double calculate_similarity(std::string seq1, std::string seq2); 

std::vector<int> create_matrix(std::string seq1, std::string seq2, 
        int rows, int cols, size_t length);

void nw_seq_to_seq( std::string& seq1, std::string& seq2, std::string& aSeq1, 
        std::string& aSeq2, std::vector<int>& M, int rows, int cols);

void calc_distances(int numSeqs, std::vector<Sequence>& seqs);

void nw_on_group( std::string& seq1, std::string& seq2, std::string& aSeq1, 
        std::string& aSeq2, std::vector<int>& M, int rows, int cols, 
        std::vector<Sequence>& group1, std::vector<Sequence>& group2);

void setup_group_alignment(std::vector<Sequence>& group1, 
        std::vector<Sequence>& group2, int g1Idx, int g2Idx); 

void choose_seq_group_align(std::vector<Sequence>& group1, 
        std::vector<Sequence>& group2); 


void align_clusters(std::vector<Sequence>& cToMerge1, 
        std::vector<Sequence>& cToMerge2); 




#endif 
