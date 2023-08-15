
#include "msa.h"
#include "matrix.h"

void calc_distances(std::vector<double>& distances, int matDims, 
        std::unique_ptr<Sequences>& seqs) {

    for (int i = 0; i < matDims; ++i) {
        for (int j = 0; j < matDims; ++j) {
            if ( i != j) {

                double dist = perform_alignment(seqs->seqs[i].seq, 
                        seqs->seqs[j].seq);

                distances[i * matDims + j] = dist;              
                //add distances to seqs
                seqs->seqs[i].distances.push_back(dist);
            } else {
                seqs->seqs[i].distances.push_back(0.0);
            }
        }
    }
}

int max_index(std::vector<double>& distances, int matDims) {
    /* Search an array and find the index of the largest element*/
    int max = -1; 
    for (int i = 0; i < (matDims * matDims); ++i) {
        if (distances[i] > max) {
            max = i; 
        }
    }

    return max; 
}


std::vector<int> create_matrix(std::string seq1, std::string seq2,
        int rows, int cols, size_t length) {
    /* create a matrix which will be traced backwards through to 
     * find the optimal sequence path.*/

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

    return M;
}


void align_seqs( std::string& seq1, std::string& seq2, std::string& aSeq1, 
        std::string& aSeq2, std::vector<int>& M, int rows, int cols) {

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

    std::vector<int> M = create_matrix(seq1, seq2, rows, cols, length);

    std::string aSeq1;
    std::string aSeq2;

    align_seqs(seq1, seq2, aSeq1, aSeq2, M, rows, cols); 

    return calculate_similarity(aSeq1, aSeq2);
}

double calculate_similarity(std::string seq1, std::string seq2) {

    int match; 
    int seqLen = seq1.length();

    for (int i = 0; i < seqLen; ++i) {
        if (seq1[i] != '-' && seq2[i] != '-' && seq1[i] == seq2[i]) {
            match++;
        }
    }

    return (double) match/seqLen;
}

void print_matrix(std::vector<double>& M, int dims) {

    for (int i = 0; i < dims; ++i) {
        for (int j = 0; j < dims; ++j) {
            std::cout << M[i * dims + j] << " ";
        }
        std::cout << std::endl; 
    }
}


void print_lower_diagnonal(std::vector<double>& lowerD, int matDims) { 

    for (int i = 1; i <= matDims; ++i) { for (int j = 1; j <= matDims; ++j) {
            if (i >= j) { //only inspect lower half of matrix 
                std::cout << lowerD[i * (i - 1) / 2 + (j - 1)] << " ";
            } else {
                std::cout << std::endl; 
                break;
            }
        }
    }
}

