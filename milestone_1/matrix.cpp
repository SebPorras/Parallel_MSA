
#include "msa.h"
#include "matrix.h"
#include <vector>

void calc_distances(int numSeqs, std::vector<Sequence>& seqs) {

    for (int i = 0; i < numSeqs; ++i) {
        for (int j = 0; j < numSeqs; ++j) {
            if ( i != j) {
                double dist = run_pairwise_alignment(seqs[i], seqs[j], false);
                //add distances to seqs
                seqs[i].distances.push_back(dist);
            } else {
                seqs[i].distances.push_back(0.0);
            }
        }
    }
}

double run_pairwise_alignment(Sequence& seq1, Sequence& seq2, bool modify) {
    /*
     * aligns two seqs and then returns their
     */

    std::string bases1 = seq1.seq;
    std::string bases2 = seq2.seq;

    //each row or column is seq length plus space for gap scores
    const int rows = bases1.length() + 1;
    const int cols = bases2.length() + 1;
    size_t length = rows * cols;

    std::vector<int> M = create_matrix(bases1, bases2, rows, cols, length);
    std::string aSeq1;
    std::string aSeq2;

    nw_seq_to_seq(bases1, bases2, aSeq1, aSeq2, M, rows, cols); 

    if (modify) { //change the actual sequences
        seq1.seq = aSeq1; 
        seq2.seq = aSeq2; 
    }

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


void nw_seq_to_seq( std::string& seq1, std::string& seq2, std::string& aSeq1, 
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

void align_clusters(std::vector<Sequence>& cToMerge1, 
        std::vector<Sequence>& cToMerge2) {

    if (cToMerge1.size() == 1 && cToMerge2.size() == 1) {
        //do a normal pairwise alignment 
        run_pairwise_alignment(cToMerge1[0], cToMerge2[0], true); 
    } else {
        choose_seq_group_align(cToMerge1, cToMerge2);
    }
}

void choose_seq_group_align(std::vector<Sequence>& group1, 
        std::vector<Sequence>& group2) {

    int g1Idx; 
    int g2Idx; 

    double mostSimiar = -1; 
    for (int i = 0; i < group1.size(); ++i) {
        for (int j = (i + 1); j < group2.size(); ++j) {

            double dist = run_pairwise_alignment(group1[i], group2[j], false); 
            if (dist > mostSimiar) {
                mostSimiar = dist; 
                g1Idx = i; 
                g2Idx = j; 
            } 
        } 
    } 

    setup_group_alignment(group1, group2, g1Idx, g2Idx);
} 

void setup_group_alignment(std::vector<Sequence>& group1, 
        std::vector<Sequence>& group2, int g1Idx, int g2Idx) {

    std::string bases1 = group1[g1Idx].seq;
    std::string bases2 = group2[g2Idx].seq;

    //each row or column is seq length plus space for gap scores
    const int rows = bases1.length() + 1;
    const int cols = bases2.length() + 1;
    size_t length = rows * cols;

    std::vector<int> M = create_matrix(bases1, bases2, rows, cols, length);

    std::string aSeq1;
    std::string aSeq2;

    nw_on_group(bases1, bases2, aSeq1, aSeq2, M, 
            rows, cols, group1, group2); 

}

void nw_on_group( std::string& seq1, std::string& seq2, std::string& aSeq1, 
        std::string& aSeq2, std::vector<int>& M, int rows, int cols, 
        std::vector<Sequence>& group1, std::vector<Sequence>& group2) {

    int I = rows - 1;  
    int J = cols - 1;   

    const int g1Size = group1.size(); 
    const int g2Size = group2.size(); 

    std::vector<std::string> g1Strs(group1.size()); 
    std::vector<std::string> g2Strs(group2.size()); 

    //same as before except now we apply changes to the whole aligned cluster
    while (I > 0 || J > 0) {

        //check left  
        if (J > 0 && M[I * cols + J] == (M[I * cols + (J - 1)] + GAP)) {

            for (int k = 0; k < g1Size; ++k) {
                //add to the front of each string 
                g1Strs[k].insert(g1Strs[k].begin(), '-');
            } 

            for (int k = 0; k < g2Size; ++k) {
                g2Strs[k].insert(g2Strs[k].begin(), group2[k].seq[J - 1]);
            } 

            J -= 1; 

            //check up  
        } if (I > 0 && M[I * cols + J] == (M[(I - 1) * cols + J] + GAP)) {

            for (int k = 0; k < g1Size; ++k) {
                g1Strs[k].insert(g1Strs[k].begin(), group1[k].seq[I - 1]);
            } 

            for (int k = 0; k < g2Size; ++k) {
                g2Strs[k].insert(g2Strs[k].begin(), '-');
            } 

            I -= 1; 

        //move diagonally 
        } else {
            for (int k = 0; k < g1Size; ++k) {
                g1Strs[k].insert(g1Strs[k].begin(), group1[k].seq[I - 1]);

            } 

            for (int k = 0; k < g2Size; ++k) {
                g2Strs[k].insert(g2Strs[k].begin(), group2[k].seq[J - 1]);
            } 

            I -= 1; 
            J -= 1; 
        }
    } 

    //update all seqs with new alignments 
    for (int k = 0; k < g1Size; ++k) {
        group1[k].seq = g1Strs[k]; 
    } 

    for (int k = 0; k < g2Size; ++k) {
        group2[k].seq = g2Strs[k]; 
    } 
}

