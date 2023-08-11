/**
 * Sebastian Porras 
 * COSC3500 
 *
 */

#include "msa.h"
#include <cstddef>
#include <cstring>
#include <iostream>
#include <iterator>

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
    int* distances = (int*)malloc(sizeof(int) * (matDims * (matDims + 1) / 2)); 

    /*
    for (int i = 0; i < matDims; ++i) {
        for (int j = 0; j < matDims; ++j) {
            if (i >= j) { //only inspect lower half of matrix 
                calc_dist(distances, i, j, seqs); 
            }
        }
    }
    */

    std::string test = "gc"; 
    std::string bar = "ga"; 

    perform_alignment(test, bar); 

    free(distances);

    return 0;
}


void calc_dist(int* distances, int i, int j, std::unique_ptr<Sequences>& seqs) {

    //dynamicAlign()
    //dists[i * (i - 1) / 2 + (j - 1) = scoreAlignment(seqI, seqJ);
}


void perform_alignment(std::string seq1, std::string seq2) {

    //each row or column is seq length plus space for gap scores
    const int rows = seq1.length() + 1;
    const int cols = seq2.length() + 1;
    size_t length = rows * cols;

    int* M = (int*)malloc(sizeof(int) * length); 
    std::memset(M, 0, sizeof(int) * length);

    int scorePenalty = 0; 
    for (int i = 0; i < cols; ++i) {
        M[i] = scorePenalty;
        scorePenalty += PENALTY; 
    }

    scorePenalty = GAP; //reset the penalty 
    for (int i = 1; i < rows; ++i) {
        
        //assign the penalty to the first column 
        M[i * rows] = scorePenalty;
        scorePenalty += PENALTY; 

        for (int j = 1; j < cols; ++j) {

            //offset seqs by one due to extra row and col 
            int diagonal = M[(i - 1) * rows + (j - 1)] 
                + (seq1[j - 1] == seq2[i - 1] ? MATCH : MISMATCH);

            int left = M[i * rows + (j - 1)] + GAP;
            int right = M[(i - 1) * rows + j] + GAP;

            //choose the best score out of our 3 directions 
            M[i * rows + j] = std::max(diagonal, std::max(left, right)); 
        }
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << M[i * rows + j] << " ";
        }
        std::cout << std::endl; 
    }


    free(M);
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

