/**
 * Sebastian Porras
 * COSC3500
 *
 */

#include "msaOpt.h"
#include "matrixOpt.h"
#include <iterator>
using namespace std;
    


int main(int argc, char **argv){

    if (argc == 1)
    {
        cout << "Provide a fasta file" << endl;
        return CLI_ERROR;
    }
    vector<Sequence> seqs = read_fasta_file(argv[FILENAME]);

    auto StartTimeRef = std::chrono::high_resolution_clock::now();

    string aOrder= "ARNDCQEGHILKMFPSTWYV";
    vector<int> subMatrix(MATRIX_SIZE, 0); 

    for (int i = 0; i < NUM_LETTERS; ++i) {
        for (int j = 0; j < NUM_LETTERS; ++j) {

            subMatrix[((int)aOrder[i] + ASCII_OFFSET) * ROW_LEN + ((int)aOrder[j] + ASCII_OFFSET)] = blosum[i][j]; 
        }
    }
  
    vector<float> distanceMatrix = calc_distances(seqs.size(), seqs, subMatrix);

    // create clusters for UPGMA
    vector<vector<Sequence>> clusters;

    for (int i = 0; i < (int) seqs.size(); ++i) {
        std::vector<Sequence> singleCluster(1, seqs[i]);
        clusters.push_back(singleCluster);
    }

    UPGMA(clusters, distanceMatrix, subMatrix);

    auto FinishTimeRef = std::chrono::high_resolution_clock::now();
    float TotalTimeRef = std::chrono::duration_cast<std::chrono::nanoseconds>(FinishTimeRef - StartTimeRef).count();
    float time = 1e-9 * TotalTimeRef;
    
    std::cout << "seconds: " << std::fixed << time << 
        std::setprecision(9) << "\n"; 
    std::cout << argv[FILENAME] << "\n"; 
    print_seqs(clusters); 

    return 0; 
}

void UPGMA(std::vector<std::vector<Sequence>> &clusters, 
        vector<float>& distanceMatrix, vector<int>& subMatrix) {

    int numClusters = clusters.size(); // iterate until there is 1 cluster
    const int numSeqs = numClusters; //track how many points we can compare

    while (numClusters > 1) {
        std::vector<Sequence> cToMerge1;
        int idxC1 = 0;

        std::vector<Sequence> cToMerge2;
        int idxC2 = 0;

        float mostSimilar = DBL_MAX;

        // find the two clusters with the
        for (int i = 0; i < numClusters; ++i) {
            for (int j = (i + 1); j < numClusters; ++j) {
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
        std::vector<Sequence> newCluster;

        for (int i = 0; i < (int) cToMerge1.size(); ++i) {
            newCluster.push_back(cToMerge1[i]);
        }

        for (int i = 0; i < (int) cToMerge2.size(); ++i) {
            newCluster.push_back(cToMerge2[i]);
        }

        clusters.push_back(newCluster);
        numClusters -= 1;
    }
}

/*
 * Find the mean difference bewteen two clusters using UPGMA.
 *
 * The difference between two clusters is defined as the
 * average difference between each pair of points to every
 * other point.https://en.wikipedia.org/wiki/UPGMA
 */
float mean_difference(std::vector<Sequence> &c1, std::vector<Sequence> &c2,
        const int numSeqs, vector<float> distanceMatrix){

    float mean = 0.0;
    const int c1Size = c1.size(); // record the size of each cluster
    const int c2Size = c2.size();
    // all sequences will have the same length for their distance array

    // take each sequence in the cluster and add up the differences
    for (int i = 0; i < c1Size; ++i)
    {
        Sequence seq1 = c1[i]; //remove loop invariants 
        int seq1Index = seq1.index;
                               //
        for (int j = 0; j < c2Size; ++j) {
            Sequence seq2 = c2[j];
            float dist = 0.0; 
            int seq2Index = seq2.index;

            for (int k = 0; k < numSeqs; ++k) {
                float delta = distanceMatrix[seq1Index * numSeqs + k] 
                    - distanceMatrix[seq2Index * numSeqs + k];                
                dist += delta * delta;
            }

            mean += std::sqrt(dist); // add the root to the mean
        }
    }

    return mean / (c1Size * c2Size);
}

/*
 * Take a fasta file and load contents into Sequence structs.
 * Will return a vector containing all the sequences. Will
 * exit with a FILE_ERROR if the file is not valid.
 */
std::vector<Sequence> read_fasta_file(std::string fileName) {
    
    std::ifstream file(fileName); // open file

    if (!file.is_open()) {
        std::cout << "File could not be opened" << std::endl;
        exit(FILE_ERROR);
    }

    std::vector<Sequence> seqs; //store sequences here 
    std::string line;
    std::string currentSeq;
    std::string currentId;
    int seqCount = 0;
    
    Sequence newSeq;
    while (std::getline(file, line))
    {

        if (line[0] == '>')
        {
            if (!currentSeq.empty()) { // save our seq        
                newSeq.seq = currentSeq; 
                newSeq.id = currentId;
                newSeq.index = seqCount;
                seqs.push_back(newSeq);
                seqCount++;

                currentSeq.clear(); // start next seq fresh
            }

            currentId = line;
        } else { // all other lines are sequences
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
 * Prints out all sequence IDs and the actual seqs
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