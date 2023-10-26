#include "msa.h"
#include "matrix.h"
using namespace std; 

//Maps the score of aligning two letters to each other 
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
void calc_distances(int numSeqs, vector<Sequence>& seqs,
                             vector<int>& subMatrix, vector<float>& distanceMatrix) {


    int worldSize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int NPerRank = int(float(numSeqs) / float(worldSize)); 
    int myFirstN = rank * NPerRank; //where in the matrix we'll start work 
    int myLastN = (rank + 1) * NPerRank; 
    
    #pragma omp parallel for collapse(2)
    for (int i = myFirstN; i < myLastN; ++i) {
        for (int j = 0; j < numSeqs; ++j) {
            //don't calculate similarity on the main diagonal 
            float dist = 0;
            if (i != j) {
                //this will return the similarity score 
                dist = run_pairwise_alignment(seqs[i], seqs[j], 
                                              false, subMatrix);
            } 

            distanceMatrix[i * numSeqs + j] = dist;
        }
    }

    MPI_Gather(&distanceMatrix[myFirstN * numSeqs], NPerRank * numSeqs, MPI_FLOAT, 
               distanceMatrix.data(), NPerRank * numSeqs, MPI_FLOAT, 0, MPI_COMM_WORLD);

    //send complete matrix back to other processes
    //MPI_Bcast(distanceMatrix.data(), numSeqs * numSeqs, MPI_FLOAT, 0, MPI_COMM_WORLD);

}

/*
 * run_pairwise_alignment
 * _______________________
 * Aligns two seqs and then returns the pairwise similarity 
 * between the two sequences. If modify is true, the sequences 
 * will be changed to their aligned version. 
 * 
 * seq1 (Sequence&): pointer to first sequence 
 * seq2 (Sequence&): pointer to second sequence 
 * modify (bool): change seq1 and seq2 to aligned version 
 * subMatrix (vector<int>&): encodes scores for aligning characters 
 *                           to one another 
 * 
 * Return (float): 
 * The similairty between the two aligned sequences
 */
float run_pairwise_alignment(Sequence& seq1, Sequence& seq2, bool modify,
                             vector<int>& subMatrix) {

    string bases1 = seq1.seq; //grab each sequence 
    string bases2 = seq2.seq;

    //each row or column is seq length plus space for gap scores
    const int rows = bases1.length() + 1;
    const int cols = bases2.length() + 1;
    size_t length = rows * cols;

    //creates the matrix which will be traced back through to find alignment 
    vector<int> M = create_matrix(bases1, bases2, rows, 
                                       cols, length, subMatrix);

    string aSeq1; //the aligned sequences will be saved here 
    string aSeq2;

    //run the NW algorithm 
    nw_seq_to_seq(bases1, bases2, aSeq1, aSeq2, M, rows, cols); 

    if (modify) { //change the actual sequences to aligned versions 
        seq1.seq = aSeq1; 
        seq2.seq = aSeq2; 
    }

    //the similarity between the two sequences 
    return calculate_similarity(aSeq1, aSeq2);
}

/*
 * nw_seq_to_seq
*  _____________
 * Walk backwards through the path matrix, M, and add gaps and chars 
 * depending on the path. Will append to new strings that are created. 
 *
 * seq1(string&): 
 * seq2 (string&):
 * aSeq1 (string&): 
 * aSeq2 (string&): , int rows, int cols
 * rows(int): len of seq A 
 * cols(int): len of seq B
 * M (vector<int>&): the path matrix 
 * 
 * ref: https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
 * 
 */
void nw_seq_to_seq(string& seq1, string& seq2, string& aSeq1, 
                   string& aSeq2, vector<int>& M, int rows, int cols) {

    int I = rows - 1;  //update the length so that you don't include gap cols
    int J = cols - 1;   

    //while you're not at the front of both sequences 
    while (I > 0 && J > 0) {

        //check if current cell matches left cell + gap score 
        if (M[I * cols + J] == (M[I * cols + (J - 1)] + GAP)) {
            
            //introduce a gap character for sequence 1 
            aSeq1 = '-' + aSeq1;
            //align the previous letter to the gap 
            aSeq2 = seq2[J - 1] + aSeq2; 
            J -= 1; 

        //check if current cell matches cell above + gap score 
        } else if (M[I * cols + J] == (M[(I - 1) * cols + J] + GAP)) {

            //align the next previous letter to the gap 
            aSeq1 = seq1[I - 1] + aSeq1;

            //introduce a gap character for sequence 2 
            aSeq2 = '-' + aSeq2;
            I -= 1; 

        //otherwise you know the best movement is to align 
        } else {
            //align the previous letters for both sequences to each other 
            aSeq1 = seq1[I - 1] + aSeq1;
            aSeq2 = seq2[J - 1] + aSeq2; 
            I -= 1; 
            J -= 1; 
        }
    } 
    
    //if one of your sequences still has characters, add gaps until you finish 
    while (I > 0) {
        aSeq1 = seq1[I - 1] + aSeq1; 
        aSeq2 = '-' + aSeq2;
        I -= 1; 
    }

    //same here
    while (J > 0) {
        aSeq1 = '-' + aSeq1;
        aSeq2 = seq2[J - 1] + aSeq2; 
        J -= 1; 
    }
}

/*
 * Calculate the pairwise similarity. It is simply 
 *  num of matching bases / seq len. 
 * 
 * seq1 (string): The first string 
 * seq2 (string): The second string 
 * 
 * Return (float):
 * The similarity score
 */
float calculate_similarity(string seq1, string seq2) {

    int match = 0; 
    int seqLen = seq1.length();

    #pragma omp parallel for reduction(+:match)
    for (int i = 0; i < seqLen; ++i) {

        //gap characters don't count so ignore them
        if (seq1[i] != '-' && seq2[i] != '-' 
            && seq1[i] == seq2[i]) {
            match++;
        }
    }

    return (float) match/seqLen;
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
 */
vector<int> create_matrix(string& seq1, string& seq2,
        const int rows, const int cols, const size_t length,
         vector<int>& subMatrix) {

    vector<int> M(length, 0); 

    #pragma omp parallel
    { 
    
    #pragma omp for schedule(static) 
    for (int i = 0; i < cols; ++i) {
        M[i] = i * GAP;
    }

    
    #pragma omp for schedule(static) 
    for (int i = 0; i < rows; ++i) {
        //assign the penalty to the first column 
        M[i * cols] = i * GAP; //avoid jumping through memory 
    }

  
    for (int I = 0; I < rows + cols - 1; I++) {
        #pragma omp for schedule(static)
        for (int J = max(0, I - rows + 1); J < min(cols, I + 1); J++) {

            int waveRow = I - J; 
            int waveCol = J;
  
            if (waveRow > 0 && waveCol > 0) {

                int diagonal = M[(waveRow - 1) * cols + (waveCol - 1)];

                //'-' chars have a score of 0, otherwise get subsitution score 
                if (seq1[waveRow - 1] != '-' &&  seq2[waveCol - 1] != '-') {
                    diagonal += subMatrix[((int)seq1[waveRow - 1] 
                    + ASCII_OFFSET) * ROW_LEN + ((int)seq2[waveCol - 1] + ASCII_OFFSET)]; 
                }

                int left = M[waveRow * cols + (waveCol - 1)] + GAP;
                int right = M[(waveRow - 1) * cols + waveCol] + GAP;
                //choose the best score out of our 3 directions 
                M[waveRow * cols + waveCol] = max(diagonal, max(left, right)); 
            }
        }
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
void align_clusters(vector<Sequence>& cToMerge1, 
        vector<Sequence>& cToMerge2, vector<int>& subMatrix) {

    if ((int) cToMerge1.size() == 1 && (int) cToMerge2.size() == 1) {
        //do a normal pairwise alignment but also modify seqs 
        run_pairwise_alignment(cToMerge1[0], cToMerge2[0], true, subMatrix); 
    } else {
        //otherwise choose which sequences to align the clustrs with 
        choose_seq_group_align(cToMerge1, cToMerge2, subMatrix);
    }
}

/*
 *choose_seq_group_align
 * ________________________
 *
 * Take two clusters, do a pairwise alignment between each sequence in 
 * the cluster and choose the most similar pair of sequences to align on.
 * Once it find sthe two most similar sequence, begin aligning the two clusters 
 *
 * group1 (vector<Sequence>&): The first cluster to check 
 * group2 (vector<Sequence>&): The second cluster being compared
 * subMatrix (vector<int>&): The matrix of alignment scores to allow alignment 
 * 
 */
void choose_seq_group_align(vector<Sequence>& group1, 
        vector<Sequence>& group2, vector<int>& subMatrix) {

    int g1Idx = 0; //allows the best sequences to be grabbed later
    int g2Idx = 0; 
    float globalClosest = -1; //similarity cannot be negative 

    #pragma omp parallel 
    {
        int localG1Idx = 0; //allows the best sequences to be grabbed later
        int localG2Idx = 0; 
        float localMostSimiar = -1; //similarity cannot be negative 

        #pragma omp for
        for (int i = 0; i < (int) group1.size(); ++i) {
            for (int j = 0; j < (int) group2.size(); ++j) {

                //calculate two pairwise alignments along and compare 
                float similarity = run_pairwise_alignment(group1[i], group2[j], 
                                                false, subMatrix); 
                

                if (similarity > localMostSimiar) {
                    localMostSimiar = similarity;
                    localG1Idx = i;
                    localG2Idx = j; 
                }

            }
        }

        #pragma omp critical
        {
            if (localMostSimiar > globalClosest) {
                
                globalClosest = localMostSimiar;
                g1Idx = localG1Idx; 
                g2Idx = localG2Idx; 
            }
        }
    }
    
    //can begin aliginng with our best two sequences from each cluster
    setup_group_alignment(group1, group2, g1Idx, g2Idx, subMatrix);
} 

/**
 * setup_group_alignment
 * _____________________
 * 
 * Same as a usual alignment, it creates the path 
 * matrix based on the alignment between 2 sequences 
 * that are most similar in these two clusters. However, 
 * it kicks off the NW algorithm on a group which updates 
 * sequences based off the alignment of the 2 sequences. 
 * 
 * group1 (vector<Sequence>&): The first cluster to check 
 * group2 (vector<Sequence>&): The second cluster being compared
 * subMatrix (vector<int>&): The matrix of alignment scores to allow alignment
 * g1SeqIdx (int): The index in the cluster vector of the seq to align on 
 * g2SeqIdx (int): The index in the cluster vector of the other seq to align on 
 * 
*/
void setup_group_alignment(vector<Sequence>& group1, 
        vector<Sequence>& group2, int g1SeqIdx, int g2SeqIdx, 
        vector<int>& subMatrix) {

    //each row or column is seq length plus space for gap scores
    const int rows = group1[g1SeqIdx].seq.length() + 1;
    const int cols = group2[g2SeqIdx].seq.length() + 1;
    const size_t length = rows * cols;

    //create the path matrix 
    vector<int> M = create_matrix(group1[g1SeqIdx].seq, group2[g2SeqIdx].seq, 
            rows, cols, length, subMatrix);

    nw_on_group(M, rows, cols, group1, group2); 
}

/**
 * nw_on_group
 * ______________
 * 
 * Performs the traceback to find optimal alignment just like nw_seq_to_seq(),
 * however, every sequence in the cluster gets gaps or alignments introduced 
 * at the same indices to preserve the previous alignments that have occured. 
 * 
 * M (vector<int>): The path matrix to trace backwards through 
 * rows (int): the length of one sequence + gap row
 * cols (int): the length of other sequence + gap columns
 * group1 (vector<Sequence>&): The first cluster to check 
 * group2 (vector<Sequence>&): The second cluster being compared
 * 
*/
void nw_on_group(vector<int>& M, int rows, int cols, 
        vector<Sequence>& group1, vector<Sequence>& group2) {

    int I = rows - 1;  //correct to make these our sequence lengths 
    int J = cols - 1;   
    const int g1Size = group1.size(); 
    const int g2Size = group2.size(); 

    vector<string> g1Strs(g1Size, ""); //these will hold alignments
    vector<string> g2Strs(g2Size, "");

    //same as before except now we apply changes to the whole aligned cluster
    while (I > 0 || J > 0) {
        //check left  
        if (J > 0 && M[I * cols + J] == (M[I * cols + (J - 1)] + GAP)) {
            
            int k; 
            for (k = 0; k < g1Size - 3; k += 4) {
                //add to the front of each string in the cluster  
                g1Strs[k] = '-' + g1Strs[k]; 
                g1Strs[k + 1] = '-' + g1Strs[k + 1]; 
                g1Strs[k + 2] = '-' + g1Strs[k + 2]; 
                g1Strs[k + 3] = '-' + g1Strs[k + 3]; 
            } 

            for (; k < g1Size; ++k) {
                g1Strs[k] = '-' + g1Strs[k]; 
            }

            int z;
            for (z = 0; z < g2Size - 3; z += 4) {
                g2Strs[z] = group2[z].seq[J - 1] + g2Strs[z];
                g2Strs[z + 1] = group2[z + 1].seq[J - 1] + g2Strs[z + 1];
                g2Strs[z + 2] = group2[z + 2].seq[J - 1] + g2Strs[z + 2];
                g2Strs[z + 3] = group2[z + 3].seq[J - 1] + g2Strs[z + 3];
            } 

            for (; z < g2Size; ++z) {
                g2Strs[z] = group2[z].seq[J - 1] + g2Strs[z];
            } 

            J -= 1; 

        //check up  
        } else if (I > 0 && M[I * cols + J] == (M[(I - 1) * cols + J] + GAP)) {
            
            int k; 
            for (k = 0; k < g1Size - 3; k += 4) {
                g1Strs[k] = group1[k].seq[I - 1] + g1Strs[k];
                g1Strs[k + 1] = group1[k + 1].seq[I - 1] + g1Strs[k + 1];
                g1Strs[k + 2] = group1[k + 2].seq[I - 1] + g1Strs[k + 2];
                g1Strs[k + 3] = group1[k + 3].seq[I - 1] + g1Strs[k + 3];
            }
            
            for (; k < g1Size; ++k) {
                g1Strs[k] = group1[k].seq[I - 1] + g1Strs[k];
            } 

            int z; 
            for (z = 0; z < g2Size - 3; z += 4) {
                g2Strs[z] = '-' + g2Strs[z]; 
                g2Strs[z + 1] = '-' + g2Strs[z + 1]; 
                g2Strs[z + 2] = '-' + g2Strs[z + 2]; 
                g2Strs[z + 3] = '-' + g2Strs[z + 3]; 
            } 

            for (; z < g2Size; ++z) {
                g2Strs[z] = '-' + g2Strs[z]; 
            } 

            I -= 1; 

        //move diagonally 
        } else {
            
            int k;
            for (k = 0; k < g1Size - 3; k += 4) {
                g1Strs[k] = group1[k].seq[I - 1] + g1Strs[k];
                g1Strs[k+ 1] = group1[k + 1].seq[I - 1] + g1Strs[k + 1];
                g1Strs[k+ 2] = group1[k + 2].seq[I - 1] + g1Strs[k + 2];
                g1Strs[k+ 3] = group1[k + 3].seq[I - 1] + g1Strs[k + 3];
            } 

            for (; k < g1Size; ++k) {
                g1Strs[k] = group1[k].seq[I - 1] + g1Strs[k];
            } 

            int z; 
            for (z = 0; z < g2Size - 3; z += 4) {
                g2Strs[z] = group2[z].seq[J - 1] + g2Strs[z];
                g2Strs[z + 1] = group2[z + 1].seq[J - 1] + g2Strs[z + 1];
                g2Strs[z + 2] = group2[z + 2].seq[J - 1] + g2Strs[z + 2];
                g2Strs[z + 3] = group2[z + 3].seq[J - 1] + g2Strs[z + 3];
            } 

            for (; z < g2Size; ++z) {
                g2Strs[z] = group2[z].seq[J - 1] + g2Strs[z];
            } 

            I -= 1; 
            J -= 1; 
        }
    } 

    //update all seqs with new alignments 
    int k; 
    for (k = 0; k < g1Size - 3; k += 4) {
        group1[k].seq = g1Strs[k];
        group1[k + 1].seq = g1Strs[k + 1];
        group1[k + 2].seq = g1Strs[k + 2];
        group1[k + 3].seq = g1Strs[k + 3];
    }

    for (; k < g1Size; ++k) {
        group1[k].seq = g1Strs[k]; 
    } 
    
    int z; 
    for (z = 0; z < g2Size - 3; z += 4) {
        group2[z].seq = g2Strs[z];
        group2[z + 1].seq = g2Strs[z + 1];
        group2[z + 2].seq = g2Strs[z + 2];
        group2[z + 3].seq = g2Strs[z + 3];
    } 

    for (; z < g2Size; ++z) {
        group2[z].seq = g2Strs[z]; 
    } 
}