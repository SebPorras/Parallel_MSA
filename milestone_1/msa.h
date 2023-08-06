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

struct Sequences {
    int numSeqs;
    std::vector<std::string> seqs;
    std::unordered_map<std::string, int> ids; 
};


void read_fasta_file(std::string fileName, std::unique_ptr<Sequences>& seqs); 

#endif
