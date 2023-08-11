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

    std::string test = "TAGC"; 
    std::string bar = "TAC"; 

    perform_alignment(bar, test); 

    free(distances);

    return 0;
}


void calc_dist(int* distances, int i, int j, std::unique_ptr<Sequences>& seqs) {

    //dynamicAlign()
    //dists[i * (i - 1) / 2 + (j - 1) = scoreAlignment(I, J);
}


void perform_alignment(std::string seq1, std::string seq2) {

    //each row or column is seq length plus space for gap scores
    const int rows = seq1.length() + 1;
    const int cols = seq2.length() + 1;
    size_t length = rows * cols;

    int* M = (int*)malloc(sizeof(int) * length); 
    std::memset(M, 0, sizeof(int) * length);


    int scorePenalty = GAP; 
    for (int i = 1; i < cols; ++i) {
        M[i] = scorePenalty;
        scorePenalty += GAP; 
    }

    scorePenalty = GAP; //reset the penalty 
    for (int i = 1; i < rows; ++i) {
        
        //assign the penalty to the first column 
        M[i * cols] = scorePenalty;
        scorePenalty += GAP; 

        for (int j = 1; j < cols; ++j) {

            //offset seqs by one due to extra row and col 
            int diagonal = M[(i - 1) * cols + (j - 1)] 
                + (seq1[i - 1] == seq2[j - 1] ? MATCH : MISMATCH);

            int left = M[i * cols + (j - 1)] + GAP;
            int right = M[(i - 1) * cols + j] + GAP;

            //choose the best score out of our 3 directions 
            M[i * cols + j] = std::max(diagonal, std::max(left, right)); 
        }
    }

   // print_matrix(M, rows, cols);

    std::string aSeq1;
    std::string aSeq2;

    int I = rows - 1;  
    int J = cols - 1;   

    while (I > 0 || J > 0) {

        //check left  
        if (J > 0 && M[I * cols + J] == (M[I * cols + (J - 1)] + GAP)) {

            aSeq1 = '-' + aSeq1;
            aSeq2 = seq2[J - 1] + aSeq2; 
            J -= 1; 

        //check up  
        } if (I > 0 && M[I * cols + J] == (M[(I - 1) * cols + J] + GAP)) {

            aSeq1 = seq1[I - 1] + aSeq1; 
            aSeq2 = '-' + aSeq2;
            I -= 1; 
            
        //move diagonally 
        } else {
            aSeq1 = seq1[I -1] + aSeq1;
            aSeq2 = seq2[J -1] + aSeq2; 
            I -= 1; 
            J -= 1; 
        }

    } 

    //std::cout << aSeq1 << std::endl; 
    //std::cout << aSeq2 << std::endl; 

    free(M);
}


void print_matrix(int* M, int rows, int cols) {

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << M[i * cols + j] << " ";
        }
        std::cout << std::endl; 
    }
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

