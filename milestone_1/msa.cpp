#include <fstream>
#include <iostream>
#include <pstl/glue_execution_defs.h>
#include "fasta_reader.cpp"

int main (int argc, char** argv) { 

    if (argc == 1) {
        std::cout << "Provide a fasta file" << std::endl;
        return 1; 
    }

    Sequences seqs; //data structure to hold our sequences 
    seqs.seqs = 0;
    seqs.numSeqs = 0;

    read_fasta_file(argv[FILENAME], &seqs);//populate with sequences

    return 0;
}
