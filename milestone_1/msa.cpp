/**
 * Sebastian Porras 
 * COSC3500 
 *
 */

#include "msa.h"
#include "matrix.h"
#include <cmath>
#include <iostream>
#include <ostream>
#include <pstl/glue_execution_defs.h>
#include <vector>

int main (int argc, char** argv) { 

    if (argc == 1) {
        std::cout << "Provide a fasta file" << std::endl;
        return CLI_ERROR; 
    }

    //data structure to hold our sequences 
    std::vector<Sequence> seqs; 

    read_fasta_file(argv[FILENAME], seqs);//populate with sequences
    calc_distances(seqs.size(), seqs);

    //create clusters for UPGMA 
    std::vector<std::vector<Sequence>> clusters; 

    for (int i = 0; i < seqs.size(); ++i) {
        std::vector<Sequence> singleCluster(1, seqs[i]);
        clusters.push_back(singleCluster);
    }
   
    UPGMA(clusters); 

    for (int i = 0; i < clusters.size(); ++i) {
        for (int j = 0; j < clusters[i].size(); j++) {
            std::cout <<  clusters[i][j].id << " : " <<
                clusters[i][j].seq << std::endl; 
        }
    } 

    return 0;
}

void UPGMA(std::vector<std::vector<Sequence>>& clusters) {

    int numClusters = clusters.size(); 
    while (numClusters > 1) {

        std::vector<Sequence> cToMerge1; 
        int idxC1;

        std::vector<Sequence> cToMerge2; 
        int idxC2;

        double mostSimilar = -1;

        for (int i = 0; i < numClusters; ++i) {
            for (int j = (i + 1); j < numClusters; ++j) {
                
                double dist = mean_difference(clusters[i], clusters[j]); 

                if (dist > mostSimilar) {
                    mostSimilar = dist; 

                    cToMerge1 = clusters[i];
                    cToMerge2 = clusters[j];

                    idxC1 = i;
                    idxC2 = j;
                }
            }
        }

        align_clusters(cToMerge1, cToMerge2); 
        
        //check which idx is greater so order is not messed up when removing 
        if (idxC1 > idxC2) {
            clusters.erase(clusters.begin() + idxC1);
            clusters.erase(clusters.begin() + idxC2);
        } else {
            clusters.erase(clusters.begin() + idxC2);
            clusters.erase(clusters.begin() + idxC1);
        }

        //collapse old clusters and remove them 
        std::vector<Sequence> newCluster;

        for (int i = 0; i < cToMerge1.size(); ++i) {
            newCluster.push_back(cToMerge1[i]);
        }

        for (int i = 0; i < cToMerge2.size(); ++i) {
            newCluster.push_back(cToMerge2[i]);
        }


        clusters.push_back(newCluster);
        numClusters -= 1; 
    }
}

double mean_difference(std::vector<Sequence>& c1, std::vector<Sequence>& c2) {
    /* The difference between two clusters is defined as the 
     * average difference between each pair of points to every 
     * other point. https://en.wikipedia.org/wiki/UPGMA*/ 
    double mean = 0.0; 
    int c1Size = c1.size();
    int c2Size = c2.size();

    for (int i = 0; i < c1Size; ++i) {
        for (int j = 0; j < c2Size; ++j) {

            Sequence seq1 = c1[i];
            Sequence seq2 = c2[j];
            int numPoints = seq1.distances.size(); 

            double dist = 0.0;
            for (int k = 0; k < numPoints; ++k) {
                float delta = seq1.distances[k] - seq2.distances[k];
                dist += delta * delta; 

            mean += std::sqrt(dist); 
            }
        }
    }

    return mean / (c1Size * c2Size);
}

void read_fasta_file(std::string fileName, std::vector<Sequence>& seqs) {
    
    std::ifstream file(fileName); //open file

    if (!file.is_open()) {
        std::cout << "File could not be opened" << std::endl; 
        exit(FILE_ERROR);
    }

    std::string line; 

    std::string currentSeq;
    std::string currentId; 
    int seqCount = 0;
    Sequence newSeq; 

    while (std::getline(file, line)) {

        if (line[0] == '>') {

            if (!currentSeq.empty()) { //save our seq 
                                       //
                newSeq.seq = currentSeq;
                newSeq.id = currentId; 
                newSeq.index = seqCount; 
                seqs.push_back(newSeq);
                seqCount++;  

                currentSeq.clear(); //start next seq fresh 
            }

            currentId = line; 

        } else { //all other lines are sequences 
            currentSeq += line; 
        }
    }

    //save the last sequence 
    newSeq.seq = currentSeq;
    newSeq.id = currentId;
    newSeq.index = seqCount; 

    seqs.push_back(newSeq);

    file.close();
}

