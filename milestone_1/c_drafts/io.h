#ifndef IO_H
#define IO_H

#include <iostream>

const int MAX_SEQ_LEN 200
const int FILENAME 1

struct FastaSeqs { 
    int numSeqs; //number of seqs we're storing
    std::array seqArray; //store the actual seq
    std::array idArray; //store the seq ID 
}; //keep track of fasta seqs we read in 

//prototypes 
void add_seq(FastaSeqs allSeqs, char* id, char* seq); 
FastaSeqs read_fasta_file(FILE* fileStream);
#endif