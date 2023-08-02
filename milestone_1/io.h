#ifndef IO_H
#define IO_H

#define MAX_SEQ_LEN 200
#define FILENAME 1

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct {
    int numSeqs; //number of seqs we're storing
    char** seqArray; //store the actual seq
    char** idArray; //store the seq ID 
} FastaSeqs; //keep track of fasta seqs we read in 

//prototypes 
void add_seq(FastaSeqs* allSeqs, char* id, char* seq); 
FastaSeqs read_fasta_file(FILE* fileStream);
#endif