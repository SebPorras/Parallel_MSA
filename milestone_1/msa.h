#ifndef MSA_H
#define MSA_H

#include <string>
#include <vector> 
#include <map>
#include <memory>
#include <list> 
#include <cmath> 

const int MAX_SEQ_LEN = 200;
const int FILENAME = 1;

const int CLI_ERROR = 1;
const int FILE_ERROR = 2;

const int MATCH = 3;
const int MISMATCH = 0;
const int GAP = -3;

struct Sequence {
    std::string seq; //actual sequence 
    std::string id; //sequence name 
    int index; //where the sequence is in the matrix  
    std::vector<double> distances; //distances to other sequences 
};

double mean_difference(std::vector<Sequence>& c1, std::vector<Sequence>& c2); 
void read_fasta_file(std::string fileName, std::vector<Sequence>& seqs); 
void UPGMA(std::vector<std::vector<Sequence>>& clusters);

#endif
