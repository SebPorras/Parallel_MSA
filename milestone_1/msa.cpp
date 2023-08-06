/**
 * Sebastian Porras 
 * COSC3500 
 *
 */

#include "msa.h"

int main (int argc, char** argv) { 

    if (argc == 1) {
        std::cout << "Provide a fasta file" << std::endl;
        return CLI_ERROR; 
    }

    //data structure to hold our sequences 
    std::unique_ptr<Sequences> seqs = std::make_unique<Sequences>();
    seqs->numSeqs = 0;

    read_fasta_file(argv[FILENAME], seqs);//populate with sequences
                                          //
    int matDims = seqs->numSeqs;
    //construct a lower diagonal matrix
    int* distances = (int*) malloc(sizeof(int) * (matDims * (matDims + 1) / 2)); 

    for (int i = 0; i < matDims; ++i) {
        for (int j = 0; j < matDims; ++j) {
            if (i >= j) { //only inspect lower half of matrix 
                calc_dist(distances, i, j, seqs); 
            }
        }
    }

    free(distances);

    return 0;
}


void calc_dist(int* distances, int i, int j, std::unique_ptr<Sequences>& seqs) {

    //dynamicAlign()
    //dists[i * (i - 1) / 2 + (j - 1) = scoreAlignment(seqI, seqJ);
}

void read_fasta_file(std::string fileName, std::unique_ptr<Sequences>& seqs) {
    
    std::ifstream file(fileName); //open file

    if (!file.is_open()) {
        std::cout << "File could not be opened" << std::endl; 
        exit(FILE_ERROR);
    }

    std::string line; 
    std::string currentSeq;
    std::string currentId; 
    while (std::getline(file, line)) {

        if (line[0] == '>') {
            seqs->ids[line] = seqs->numSeqs; //map id to index 
            seqs->numSeqs++; 
            if (!currentSeq.empty()) { //save our seq 
                seqs->seqs.push_back(currentSeq);
                currentSeq.clear();
            }

        } else { //all other lines are sequences 
            currentSeq += line; 
        }
    }

    //save the last sequence 
    seqs->seqs.push_back(currentSeq);
    file.close();
}

