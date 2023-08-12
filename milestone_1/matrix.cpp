
#include "msa.h"
#include "matrix.h"

void calc_distances(std::vector<double>& distances, int matDims, 
        std::unique_ptr<Sequences>& seqs) {

    for (int i = 1; i <= matDims; ++i) {
        for (int j = 1; j <= matDims; ++j) {
            if (i >= j) { //only inspect lower half of matrix 
                distances[i * (i - 1) 
                    / 2 + (j - 1)] 
                    = perform_alignment(seqs->seqs[i - 1], seqs->seqs[j - 1]);
            }
        }
    }
}

void print_lower_diagnonal(std::vector<double>& lowerD, int matDims) {

    for (int i = 1; i <= matDims; ++i) {
        for (int j = 1; j <= matDims; ++j) {
            if (i >= j) { //only inspect lower half of matrix 
                std::cout << lowerD[i * (i - 1) / 2 + (j - 1)] << " ";
            } else {
                std::cout << std::endl; 
                break;
            }
        }
    }
}


/*
 *
 * aligns two seqs and then returns their
 */
double perform_alignment(std::string seq1, std::string seq2) {

    //each row or column is seq length plus space for gap scores
    const int rows = seq1.length() + 1;
    const int cols = seq2.length() + 1;
    size_t length = rows * cols;

    std::vector<int> M(length, 0); 

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

    return calculate_similarity(aSeq1, aSeq2);
}


double calculate_similarity(std::string seq1, std::string seq2) {

    int matches; 
    int seqLen = seq1.length();

    for (int i = 0; i < seqLen; ++i) {
        if (seq1[i] != '-' && seq2[i] != '-' && seq1[i] == seq2[i]) {
            matches++;
        }
    }

    return (double) matches/seqLen;
}



void print_matrix(int* M, int rows, int cols) {

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << M[i * cols + j] << " ";
        }
        std::cout << std::endl; 
    }
}


