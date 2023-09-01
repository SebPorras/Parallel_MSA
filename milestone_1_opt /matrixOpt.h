#ifndef MATRIXOPT_H 
#define MATRIXOPT_H

#include "msaOpt.h"
#include <ostream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstddef>
#include <iterator>
#include <unordered_map>
using namespace std; 

float run_pairwise_alignment(Sequence& seq1, Sequence& seq2, bool modify, vector<int>& subMatrix);

float calculate_similarity(std::string seq1, std::string seq2); 

std::vector<int> create_matrix(std::string& seq1, std::string& seq2, 
        const int rows, const int cols, const size_t length, vector<int>& subMatrix);

void nw_seq_to_seq( std::string& seq1, std::string& seq2, std::string& aSeq1, 
        std::string& aSeq2, std::vector<int>& M, int rows, int cols);

vector<float> calc_distances(int numSeqs, std::vector<Sequence>& seqs, vector<int>& subMatrix);

void nw_on_group(std::vector<int>& M, int rows, int cols, 
        std::vector<Sequence>& group1, std::vector<Sequence>& group2);

void setup_group_alignment(std::vector<Sequence>& group1, 
        std::vector<Sequence>& group2, int g1Idx, int g2Idx, 
        vector<int>& subMatrix); 

void choose_seq_group_align(std::vector<Sequence>& group1, 
        std::vector<Sequence>& group2, vector<int>& subMatrix); 


void align_clusters(std::vector<Sequence>& cToMerge1, 
        std::vector<Sequence>& cToMerge2, vector<int>& subMatrix); 


                

#endif 
