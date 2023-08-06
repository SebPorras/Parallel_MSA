/**
 * Sebastian Porras 
 * COSC3500 
 *
 */

#include "msa.h"

int main (int argc, char** argv) { 

    if (argc == 1) {
        std::cout << "Provide a fasta file" << std::endl;
        return 1; 
    }

    //data structure to hold our sequences 
    std::unique_ptr<Sequences> seqs = std::make_unique<Sequences>();
    seqs->numSeqs = 0;

    read_fasta_file(argv[FILENAME], seqs);//populate with sequences
    std::cout << "num seqs " << seqs->numSeqs << std::endl;
    
    for (const auto& pair : seqs->ids) {
        std::cout << pair.first << " " << pair.second << std::endl; 
        std::cout << seqs->seqs[pair.second] << std::endl; 
    }

    return 0;
}

void read_fasta_file(std::string fileName, std::unique_ptr<Sequences>& seqs) {
    
    std::ifstream file(fileName); //open file

    if (!file.is_open()) {
        std::cout << "File could not be opened" << std::endl; 
        exit(2);
    }

    std::string line; 
    std::string currentSeq;
    std::string currentId; 

    while (std::getline(file, line)) {

        if (line[0] == '>') {
            seqs->ids[line] = seqs->numSeqs; //map id to index 
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

