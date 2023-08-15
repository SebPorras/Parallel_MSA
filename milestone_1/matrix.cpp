
#include "msa.h"
#include "matrix.h"
#include <vector>

void calc_distances(int numSeqs, std::vector<Sequence>& seqs) {

    for (int i = 0; i < numSeqs; ++i) {
        for (int j = 0; j < numSeqs; ++j) {
            if ( i != j) {
                double dist = perform_alignment(seqs[i], seqs[j], false);
                //add distances to seqs
                seqs[i].distances.push_back(dist);
            } else {
                seqs[i].distances.push_back(0.0);
            }
        }
    }
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

            aSeq1 = 'X' + aSeq1;
            aSeq2 = seq2[J - 1] + aSeq2; 
            J -= 1; 

            //check up  
        } if (I > 0 && M[I * cols + J] == (M[(I - 1) * cols + J] + GAP)) {

            aSeq1 = seq1[I - 1] + aSeq1; 
            aSeq2 = 'X' + aSeq2;
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
double perform_alignment(Sequence& seq1, Sequence& seq2, bool modify) {

    std::string bases1 = seq1.seq;
    std::string bases2 = seq2.seq;

    //each row or column is seq length plus space for gap scores
    const int rows = bases1.length() + 1;
    const int cols = bases2.length() + 1;
    size_t length = rows * cols;

    std::vector<int> M = create_matrix(bases1, bases2, rows, cols, length);

    std::string aSeq1;
    std::string aSeq2;

    align_seqs(bases1, bases2, aSeq1, aSeq2, M, rows, cols); 

    if (modify) { //change the actual sequences
        seq1.seq = aSeq1; 
        seq2.seq = aSeq2; 
    }

    return calculate_similarity(aSeq1, aSeq2);
}


void align_seq_to_group(std::vector<Sequence> seq, 
        std::vector<Sequence> group) {


    return; 
} 

void align_group_to_group(std::vector<Sequence> group1, 
        std::vector<Sequence> group2) {
    

    return; 
} 


void align_clusters(std::vector<Sequence> cToMerge1, 
        std::vector<Sequence> cToMerge2) {

    if (cToMerge1.size() == 1 && cToMerge2.size() == 1) {
        //do a normal pairwise alignment 
        perform_alignment(cToMerge1[0], cToMerge2[0], true); 

    } else if (cToMerge1.size() == 1 && cToMerge2.size() > 1) {
        align_seq_to_group(cToMerge1, cToMerge2);

    } else if (cToMerge1.size() > 1 && cToMerge2.size() == 1) {
        align_seq_to_group(cToMerge2, cToMerge1);
    }

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

