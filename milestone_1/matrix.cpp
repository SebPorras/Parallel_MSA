
#include "msa.h"
#include "matrix.h"
#include <vector>

int blosum[20][20] = {
 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0,
-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3,
-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,
-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,
 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,
-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,
-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,
 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3,
-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,
-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3,
-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1,
-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,
-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1,
-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1,
-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2,
 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,
 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0,
-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3,
-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1,
 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4
 };

std::unordered_map<char, int> acids = {
    {'A', 0}, {'R', 1}, {'N', 2}, {'D', 3}, {'C', 4}, {'Q', 5},
    {'E', 6}, {'G', 7}, {'H', 8}, {'I', 9}, {'L', 10}, {'K', 11},
    {'M', 12}, {'F', 13}, {'P', 14}, {'S', 15}, {'T', 16}, {'W', 17},
    {'Y', 18}, {'V', 19}
};

/*
 * calc_distances
 * ______________
 *
 * For each sequence in the vector, runs a 
 * pairwise alignment (NW algorithm) and save the distance 
 * to a vector within the Sequence struct which saves all 
 * pairwise distances between other sequences. This will 
 * be used to cluster on when deciding the order of alignment. 
 *
 * numSeqs (int): the number of sequences to be aligned 
 * seqs (vector): the array of all Sequence structs 
 */
void calc_distances(int numSeqs, std::vector<Sequence>& seqs) {

    for (int i = 0; i < numSeqs; ++i) {
        for (int j = 0; j < numSeqs; ++j) {
            if ( i != j) {
                double dist = run_pairwise_alignment(seqs[i], seqs[j], false);
                //add distances to seqs
                seqs[i].distances.push_back(dist);
            } else {
                //ignore the diagonal 
                seqs[i].distances.push_back(0.0);
            }
        }
    }
}

/*
 * run_pairwise_alignment
 * _______________________
 * Aligns two seqs and then returns the pairwise similarity 
 * between the two sequences. If modify is true, the sequences 
 * will be changed to their aligned version. 
 */
double run_pairwise_alignment(Sequence& seq1, Sequence& seq2, bool modify) {

    std::string bases1 = seq1.seq;
    std::string bases2 = seq2.seq;

    //each row or column is seq length plus space for gap scores
    const int rows = bases1.length() + 1;
    const int cols = bases2.length() + 1;
    size_t length = rows * cols;

    //creates the matrix which will be traced back through to find alignment 
    std::vector<int> M = create_matrix(bases1, bases2, rows, cols, length);
    std::string aSeq1; //the aligned sequences 
    std::string aSeq2;

    //run the NW algorithm 
    nw_seq_to_seq(bases1, bases2, aSeq1, aSeq2, M, rows, cols); 

    if (modify) { //change the actual sequences
        seq1.seq = aSeq1; 
        seq2.seq = aSeq2; 
    }

    //the % similarity between the two sequences 
    return calculate_similarity(aSeq1, aSeq2);
}

/*
 * Walk backwards through the path matrix, M, and add gaps and chars 
 * depending on the path. Will append to new strings that are created. 
 *
 * M (vector<int>) : the path matrix 
 * rows: len of seq A 
 * cols: len of seq B 
 */
void nw_seq_to_seq(std::string& seq1, std::string& seq2, std::string& aSeq1, 
        std::string& aSeq2, std::vector<int>& M, int rows, int cols) {

    int I = rows - 1;  
    int J = cols - 1;   

    while (I > 0 && J > 0) {

        //check left  
        if (M[I * cols + J] == (M[I * cols + (J - 1)] + GAP)) {

            aSeq1 = '-' + aSeq1;
            aSeq2 = seq2[J - 1] + aSeq2; 
            J -= 1; 

            //check up  
        } else if (M[I * cols + J] == (M[(I - 1) * cols + J] + GAP)) {

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

    while (I > 0) {
        aSeq1 = seq1[I - 1] + aSeq1; 
        aSeq2 = '-' + aSeq2;
        I -= 1; 
    }

    while (J > 0) {
        aSeq1 = '-' + aSeq1;
        aSeq2 = seq2[J - 1] + aSeq2; 
        J -= 1; 
    }
}


/*
 * Calculate the pairwise similarity. It is simply 
 *  num of matching bases / seq len. 
 */
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


int get_sub_score(char i, char j) {
    return  blosum[acids[i]][acids[j]];
}

/* create a matrix which will be traced backwards through to 
 * find the optimal sequence path.
 * 
 * rows: the length of sequence 1 
 * cols: the length of sequence 2 
 * length: rows * cols 
 *
 * Return a vector with a length of rows * cols filled 
 * with scores for all possible paths through the matrix. 
 * */

std::vector<int> create_matrix(std::string& seq1, std::string& seq2,
        const int rows, const int cols, const size_t length) {

    std::vector<int> M(length, 0); 

    int scorePenalty = GAP; //top row has all gaps 
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

            //offset seqs by one due to extra row and col for gaps
            int diagonal = M[(i - 1) * cols + (j - 1)] 
            + get_sub_score(seq1[i - 1], seq2[j - 1]); 
            //    + (seq1[i - 1] == seq2[j - 1] ? MATCH : MISMATCH);

            int left = M[i * cols + (j - 1)] + GAP;
            int right = M[(i - 1) * cols + j] + GAP;

            //choose the best score out of our 3 directions 
            M[i * cols + j] = std::max(diagonal, std::max(left, right)); 
        }
    }
    return M;
}
/*
 * Take two clusters of sequences and align them. If the clusters 
 * are sinlge seqs, simply perform the NW alignment. Otherwise, every 
 * sequence in the cluster will be compared to determine the most suitable 
 * alignment before aligning on that sequence. 
 */
void align_clusters(std::vector<Sequence>& cToMerge1, 
        std::vector<Sequence>& cToMerge2) {

    if ((int) cToMerge1.size() == 1 && (int) cToMerge2.size() == 1) {
        //do a normal pairwise alignment but also modify seqs 
        run_pairwise_alignment(cToMerge1[0], cToMerge2[0], true); 
    } else {
        choose_seq_group_align(cToMerge1, cToMerge2);
    }
    
}

/*
 * Take two clusters, do a pairwise alignment between each sequence in 
 * the cluster and choose the most similar pair of sequences to align on. 
 *
 */
void choose_seq_group_align(std::vector<Sequence>& group1, 
        std::vector<Sequence>& group2) {

    int g1Idx; //allows the best sequences to be grabbed later
    int g2Idx; 

    double mostSimiar = -1; //similarity cannot be negative 
    for (int i = 0; i < (int) group1.size(); ++i) {
        for (int j = (i + 1); j < (int) group2.size(); ++j) {

            double dist = run_pairwise_alignment(group1[i], group2[j], false); 

            //update if we find a sequence that is more similar
            if (dist > mostSimiar) {
                mostSimiar = dist; 
                g1Idx = i; 
                g2Idx = j; 
            } 
        } 
    } 

    //can begin aliginng with our best two sequences from each cluster
    setup_group_alignment(group1, group2, g1Idx, g2Idx);
} 

void setup_group_alignment(std::vector<Sequence>& group1, 
        std::vector<Sequence>& group2, int g1Idx, int g2Idx) {

    //each row or column is seq length plus space for gap scores
    const int rows = group1[g1Idx].seq.length() + 1;
    const int cols = group2[g2Idx].seq.length() + 1;
    const size_t length = rows * cols;

    //create the path matrix 
    std::vector<int> M = create_matrix(group1[g1Idx].seq, group2[g2Idx].seq, 
            rows, cols, length);

    nw_on_group(M, rows, cols, group1, group2); 
}

void nw_on_group(std::vector<int>& M, int rows, int cols, 
        std::vector<Sequence>& group1, std::vector<Sequence>& group2) {


    int I = rows - 1;  
    int J = cols - 1;   
    const int g1Size = group1.size(); 
    const int g2Size = group2.size(); 

    std::vector<std::string> g1Strs(g1Size, ""); //these will hold alignments
    std::vector<std::string> g2Strs(g2Size, "");

    //same as before except now we apply changes to the whole aligned cluster
    while (I > 0 || J > 0) {
        //check left  
        if (J > 0 && M[I * cols + J] == (M[I * cols + (J - 1)] + GAP)) {

            for (int k = 0; k < g1Size; ++k) {
                //add to the front of each string in the cluster  
                g1Strs[k] = '-' + g1Strs[k]; 
            } 

            for (int k = 0; k < g2Size; ++k) {
                g2Strs[k] = group2[k].seq[J - 1] + g2Strs[k];
            } 
            J -= 1; 

            //check up  
        } else if (I > 0 && M[I * cols + J] == (M[(I - 1) * cols + J] + GAP)) {

            for (int k = 0; k < g1Size; ++k) {
                g1Strs[k] = group1[k].seq[I - 1] + g1Strs[k];
            } 

            for (int k = 0; k < g2Size; ++k) {
                g2Strs[k] = '-' + g2Strs[k]; 
            } 
            I -= 1; 

        //move diagonally 
        } else {
            for (int k = 0; k < g1Size; ++k) {
                g1Strs[k] = group1[k].seq[I - 1] + g1Strs[k];
            } 

            for (int k = 0; k < g2Size; ++k) {
                g2Strs[k] = group2[k].seq[J - 1] + g2Strs[k];
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

