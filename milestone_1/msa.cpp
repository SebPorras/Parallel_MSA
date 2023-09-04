/**
 * Sebastian Porras
 * COSC3500 - Milestone 1 
 * 
 * An implementation of a progressive sequence 
 * alignment for protein sequences. 
 */

#include "msa.h"
#include "matrix.h"
using namespace std;
    
int main(int argc, char **argv){

    if (argc == 1) {
        cout << "Provide a fasta file" << endl;
        return CLI_ERROR;
    }

    vector<Sequence> seqs = read_fasta_file(argv[FILENAME]);

    auto StartTimeRef = chrono::high_resolution_clock::now();

    //convert blosum into a direct access array 
    vector<int> subMatrix = make_sub_matrix(); 

    //calculate similarity matrix between all pairs of sequences 
    vector<float> distanceMatrix = calc_distances(seqs.size(), seqs, subMatrix);

    // Assign each sequence to its own cluster
    vector<vector<Sequence>> clusters;
    for (int i = 0; i < (int) seqs.size(); ++i) {
        vector<Sequence> singleCluster(1, seqs[i]);
        clusters.push_back(singleCluster);
    }

    //perform the clustering, aligning clusters as you go 
    UPGMA(clusters, distanceMatrix, subMatrix);

    auto FinishTimeRef = chrono::high_resolution_clock::now();
    float TotalTimeRef = chrono::duration_cast<chrono::nanoseconds>(FinishTimeRef - StartTimeRef).count();
    float time = 1e-9 * TotalTimeRef;
    
    cout << "seconds: " << fixed << time << 
        setprecision(9) << "\n"; 
    cout << argv[FILENAME] << "\n"; 

    print_seqs(clusters); 

    return 0; 
}

/**
 * Go through the blosum matrix and 
 * convert it into a direct access array 
 * based on ASCII characters of the different 
 * possible chemicals. 
 * 
 * Return (vector<int>):
 * The direct access array 
*/
vector<int> make_sub_matrix(void) {

    //the order of chemicals stored in the blosum matirx 
    string aOrder= "ARNDCQEGHILKMFPSTWYV";
    vector<int> subMatrix(MATRIX_SIZE, 0); 

    //map the substitution scores based on ASCII values 
    for (int i = 0; i < NUM_LETTERS; i += 4) {
        for (int j = 0; j < NUM_LETTERS; j += 4) {
            
            subMatrix[((int)aOrder[i] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j] + ASCII_OFFSET)] = blosum[i][j]; 
            
            subMatrix[((int)aOrder[i] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 1] + ASCII_OFFSET)] = blosum[i][j + 1]; 
            
            subMatrix[((int)aOrder[i] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 2] + ASCII_OFFSET)] = blosum[i][j + 2]; 
            
            subMatrix[((int)aOrder[i] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 3] + ASCII_OFFSET)] = blosum[i][j + 3]; 
        }

        for (int j = 0; j < NUM_LETTERS; j += 4) {
            subMatrix[((int)aOrder[i + 1] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j] + ASCII_OFFSET)] = blosum[i+ 1][j]; 

            subMatrix[((int)aOrder[i + 1] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 1] + ASCII_OFFSET)] = blosum[i + 1][j + 1]; 
            
            subMatrix[((int)aOrder[i + 1] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 2] + ASCII_OFFSET)] = blosum[i + 1][j + 2]; 
            
            subMatrix[((int)aOrder[i + 1] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 3] + ASCII_OFFSET)] = blosum[i + 1][j + 3]; 
        }

        for (int j = 0; j < NUM_LETTERS; j += 4) {
            subMatrix[((int)aOrder[i + 2] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j] + ASCII_OFFSET)] = blosum[i + 2][j]; 

            subMatrix[((int)aOrder[i + 2] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 1] + ASCII_OFFSET)] = blosum[i + 2][j + 1]; 

            subMatrix[((int)aOrder[i + 2] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 2] + ASCII_OFFSET)] = blosum[i + 2][j + 2]; 

            subMatrix[((int)aOrder[i + 2] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 3] + ASCII_OFFSET)] = blosum[i + 2][j + 3]; 
        }

        for (int j = 0; j < NUM_LETTERS; j += 4) {
            
            subMatrix[((int)aOrder[i + 3] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j] + ASCII_OFFSET)] = blosum[i + 3][j]; 
            
            subMatrix[((int)aOrder[i + 3] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 1] + ASCII_OFFSET)] = blosum[i + 3][j + 1]; 
            
            subMatrix[((int)aOrder[i + 3] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 2] + ASCII_OFFSET)] = blosum[i + 3][j + 2]; 
            
            subMatrix[((int)aOrder[i + 3] + ASCII_OFFSET) * ROW_LEN 
            + ((int)aOrder[j + 3] + ASCII_OFFSET)] = blosum[i + 3][j + 3];  
        }
    }

    return subMatrix; 
}

void UPGMA(vector<vector<Sequence>> &clusters, 
        vector<float>& distanceMatrix, vector<int>& subMatrix) {

    int numClusters = clusters.size(); // iterate until there is 1 cluster
    const int numSeqs = numClusters; //track how many points we can compare

    while (numClusters > 1) {

        vector<Sequence> cToMerge1;
        int idxC1 = 0;
        vector<Sequence> cToMerge2;
        int idxC2 = 0;
        float mostSimilar = DBL_MAX;

        // find the two closest clusters 
        for (int i = 0; i < numClusters; ++i) {
            int j; 
            for (j = i + 1; j < numClusters - 1; j += 2) {

                float dist_1 = mean_difference(clusters[i], clusters[j], 
                        numSeqs, distanceMatrix);
                
                float dist_2 = mean_difference(clusters[i], clusters[j + 1], 
                numSeqs, distanceMatrix);

                float small_chunk = (dist_1 < dist_2) ? dist_1 : dist_2;
                float small_j = (dist_1 < dist_2) ? j : j + 1; 

                if (small_chunk < mostSimilar) {
                    mostSimilar = small_chunk;

                    cToMerge1 = clusters[i];
                    cToMerge2 = clusters[small_j];

                    idxC1 = i;
                    idxC2 = small_j;
                }
            }
            
            for (; j < numClusters; ++j) {
                float dist = mean_difference(clusters[i], clusters[j], 
                        numSeqs, distanceMatrix);

                if (dist < mostSimilar) {
                    mostSimilar = dist;

                    cToMerge1 = clusters[i];
                    cToMerge2 = clusters[j];

                    idxC1 = i;
                    idxC2 = j;
                }
            }
        }

        align_clusters(cToMerge1, cToMerge2, subMatrix);

        // check which idx is greater so order is not messed up when removing
        if (idxC1 > idxC2) {
            clusters.erase(clusters.begin() + idxC1);
            clusters.erase(clusters.begin() + idxC2);
        } else {
            clusters.erase(clusters.begin() + idxC2);
            clusters.erase(clusters.begin() + idxC1);
        }

        // collapse old clusters and remove them
        vector<Sequence> newCluster;

        int i; 
        for (i = 0; i < (int) cToMerge1.size() - 3; i += 4) {
            newCluster.push_back(cToMerge1[i]);
            newCluster.push_back(cToMerge1[i + 1]);
            newCluster.push_back(cToMerge1[i + 2]);
            newCluster.push_back(cToMerge1[i + 3]);
        }
        
        for (; i < (int) cToMerge1.size(); ++i) {
            newCluster.push_back(cToMerge1[i]);
        }

        for (i = 0; i < (int) cToMerge2.size() - 3; i += 4) {
            newCluster.push_back(cToMerge2[i]);
            newCluster.push_back(cToMerge2[i + 1]);
            newCluster.push_back(cToMerge2[i + 2]);
            newCluster.push_back(cToMerge2[i + 3]);
        }

        for (; i < (int) cToMerge2.size(); ++i) {
             newCluster.push_back(cToMerge2[i]);
        }

        clusters.push_back(newCluster);
        numClusters -= 1;
    }
}

void find_closest_clusters(int numClusters, vector<vector<Sequence>> &clusters,
                           int numSeqs, vector<float>& distanceMatrix, 
                           vector<Sequence>* cToMerge1, int* idxC1, 
                           vector<Sequence>* cToMerge2, int* idxC2,
                           float* mostSimilar) {
 
    for (int i = 0; i < numClusters; ++i) {
        int j; 
        for (j = i + 1; j < numClusters - 1; j += 2) {

            float dist_1 = mean_difference(clusters[i], clusters[j], 
                    numSeqs, distanceMatrix);
            
            float dist_2 = mean_difference(clusters[i], clusters[j + 1], 
            numSeqs, distanceMatrix);

            float small_chunk = (dist_1 < dist_2) ? dist_1 : dist_2;
            float small_j = (dist_1 < dist_2) ? j : j + 1; 

            if (small_chunk < *mostSimilar) {
                *mostSimilar = small_chunk;

                *cToMerge1 = clusters[i];
                *cToMerge2 = clusters[small_j];

                *idxC1 = i;
                *idxC2 = small_j;
            }
        }
        
        for (; j < numClusters; ++j) {
            float dist = mean_difference(clusters[i], clusters[j], 
                    numSeqs, distanceMatrix);

            if (dist < *mostSimilar) {
                *mostSimilar = dist;

                *cToMerge1 = clusters[i];
                *cToMerge2 = clusters[j];

                *idxC1 = i;
                *idxC2 = j;
            }
        }
    }
}


/*
 * Find the mean difference bewteen two clusters using UPGMA.
 *
 * The difference between two clusters is defined as the
 * average difference between each pair of points to every
 * other point. https://en.wikipedia.org/wiki/UPGMA
 * 
 * c1 (vector<Sequence>&): the first cluster
 * c2 (vector<Sequence>&): the second cluster 
 * numSeqs (int): The number of sequences in the distance array 
 * distanceMatrix (vector<float>): records the distances between all sequences
 * 
 * Return (float): 
 * 
 */
float mean_difference(vector<Sequence> &c1, vector<Sequence> &c2,
        const int numSeqs, vector<float> distanceMatrix){

    float mean = 0.0; //will store average distance 
    const int c1Size = c1.size(); // record the size of each cluster
    const int c2Size = c2.size();

    // take each sequence in the cluster and add up the differences
    for (int i = 0; i < c1Size; ++i)
    {
        Sequence seq1 = c1[i]; //remove loop invariants 
        int seq1Index = seq1.index; //will be used to index into dist matrix
                              
        for (int j = 0; j < c2Size; ++j) {
            Sequence seq2 = c2[j];
            float dist = 0.0; 
            int seq2Index = seq2.index;

            //iterate through the distance matrix and compare similarity 
            // to other sequences 
            int k; 
            for (k = 0; k < numSeqs - 3; k += 4) {
                float delta = distanceMatrix[seq1Index * numSeqs + k] 
                              - distanceMatrix[seq2Index * numSeqs + k];                

                float delta2 = distanceMatrix[seq1Index * numSeqs + k + 1] 
                               - distanceMatrix[seq2Index * numSeqs + k + 1];                
             
                float delta3 = distanceMatrix[seq1Index * numSeqs + k + 2]
                               - distanceMatrix[seq2Index * numSeqs + k + 2];                
            
                float delta4 = distanceMatrix[seq1Index * numSeqs + k + 3] 
                               - distanceMatrix[seq2Index * numSeqs + k + 3]; 
                               
                dist += (delta * delta + delta2 * delta2 
                         + delta3 * delta3 + delta4 * delta4);
            }

            //add up any remaining differences 
            for (; k < numSeqs; ++k) {
                float delta = distanceMatrix[seq1Index * numSeqs + k] 
                              - distanceMatrix[seq2Index * numSeqs + k];                
                dist += delta * delta;
            }

            mean += sqrt(dist); // add the square root to the mean
        }
    }

    return mean / (c1Size * c2Size);
}

/*
 * Take a fasta file and load contents into Sequence structs.
 * Will return a vector containing all the sequences. Will
 * exit with a FILE_ERROR if the file is not valid.
 * 
 * fileName(string): name of fasta file
 */
vector<Sequence> read_fasta_file(string fileName) {
    
    ifstream file(fileName); // open file

    if (!file.is_open()) {
        cout << "File could not be opened" << endl;
        exit(FILE_ERROR);
    }

    vector<Sequence> seqs; //store sequences here 
    string line; //current line we're on 
    string currentSeq;
    string currentId;
    int seqCount = 0;
    
    Sequence newSeq;
    while (getline(file, line)) {

        //indicates new sequence 
        if (line[0] == '>') {

            // our current seq is done, save it and clear space for next one 
            if (!currentSeq.empty()) {        
                newSeq.seq = currentSeq; 
                newSeq.id = currentId;
                newSeq.index = seqCount;
                seqs.push_back(newSeq);
                seqCount++;

                currentSeq.clear(); // start next seq fresh
            }
            
            //our new ID is the current line 
            currentId = line;

        } else { // all other lines are sequence information 
            currentSeq += line;
        }
    }
    // save the last sequence
    newSeq.seq = currentSeq;
    newSeq.id = currentId;
    newSeq.index = seqCount;

    seqs.push_back(newSeq);
    file.close();

    return seqs;
}

/**
 * Prints out sequence IDs followed by actual sequence 
 * as per the FASTA format
*/
void print_seqs(vector<vector<Sequence>> clusters) {
    for (int i = 0; i < (int) clusters.size(); ++i)
    {
        for (int j = 0; j < (int) clusters[i].size(); j++)
        {
            cout << clusters[i][j].id << "\n" << 
                clusters[i][j].seq << endl;
        }
    }
}