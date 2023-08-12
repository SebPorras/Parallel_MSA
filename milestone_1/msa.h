#ifndef MSA_H
#define MSA_H

#include <string>
#include <vector> 
#include <map>
#include <memory>

const int MAX_SEQ_LEN = 200;
const int FILENAME = 1;

const int CLI_ERROR = 1;
const int FILE_ERROR = 2;


const int MATCH = 3;
const int MISMATCH = 0;
const int GAP = -10;

struct Sequences {
    int numSeqs;
    std::vector<std::string> seqs;
    std::map<std::string, int> ids; 
    std::map<int, std::string> idToName; 
};

void read_fasta_file(std::string fileName, std::unique_ptr<Sequences>& seqs); 
#endif
