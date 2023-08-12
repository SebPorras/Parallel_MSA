/**
 * Sebastian Porras 
 * COSC3500 
 *
 */

#include "msa.h"
#include "matrix.h"

int main (int argc, char** argv) { 

    if (argc == 1) {
        std::cout << "Provide a fasta file" << std::endl;
        return CLI_ERROR; 
    }

    //data structure to hold our sequences 
    std::unique_ptr<Sequences> seqs = std::make_unique<Sequences>();
    seqs->numSeqs = 0;

    read_fasta_file(argv[FILENAME], seqs);//populate with sequences
                                          
    int matDims = seqs->numSeqs;

    //construct a lower diagonal matrix
    std::vector<double> distances((matDims * (matDims + 1) / 2), 0.0); 
    calc_distances(distances, matDims, seqs);

    print_lower_diagnonal(distances, matDims);

    return 0;
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
            seqs->idToName[seqs->numSeqs] = line; //map id to index 
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

