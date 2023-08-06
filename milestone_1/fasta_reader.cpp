
#include <fstream>
#include <iostream>
#include <map>
#include <cstdio>
#include <ostream>
#include <string>
#include <cstring>

const int MAX_SEQ_LEN = 200;
const int FILENAME = 1;

struct Sequences {
    int numSeqs;
    char* seqs;
    std::map<std::string, int> ids; 
};


void read_fasta_file(std::string fileName, Sequences* seqs) {
    

    std::ifstream file(fileName); //open file

    if (!file.is_open()) {
        std::cout << "File could not be opened" << std::endl; 
        exit(2);
    }

    std::string line; 
    while (std::getline(file, line)) {
        if (line[0] == '>') {
            seqs->ids[line] = seqs->numSeqs;
        }

    }

    file.close();
}

