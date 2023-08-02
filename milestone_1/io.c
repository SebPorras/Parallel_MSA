#include "io.h"

/*
 * Read in all seqs from a fasta file 
 * and save them to the a FastaSeqs struct 
 * for reference.
 * 
 * fileStream: the fasta filename.
 * 
 * Return: A FastaSeqs struct of all seqs
*/
FastaSeqs read_fasta_file(FILE* fileStream) {

    FastaSeqs allSeqs; 
    allSeqs.numSeqs = 0;
    allSeqs.seqArray = 0;
    allSeqs.idArray = 0;

    char currentId[MAX_SEQ_LEN]; //buffer for next id
    char currentSeq[MAX_SEQ_LEN]; //read the next line which will have a sequence 
    while (fgets(currentId, MAX_SEQ_LEN, fileStream)) {


        if (currentId[0] == '>') {

            //null terminate any line that will be saved
            currentId[strlen(currentId) - 1] = '\0';

            fgets(currentSeq, MAX_SEQ_LEN, fileStream); 
            currentSeq[strlen(currentSeq) - 1] = '\0';

            add_seq(&allSeqs, currentId, currentSeq);
        }
    }

    return allSeqs;
}


/**
 * Allocates memory to store a Fasta ID and sequence 
 * and stores it within a FastaSeqs struct.
 * 
 * allSeqs: the FastaSeqs struct to store the new sequence 
 * id: the ID of the new fasta seq 
 * seq: the actual sequence 
 * 
 * Return: the FastaSeq struct with the new entry 
*/
void add_seq(FastaSeqs* allSeqs, char* id, char* seq) {

    char* idCopy;
    char* seqCopy;
    idCopy = strdup(id);
    seqCopy = strdup(seq);
            
    //add extra space in the two arrays 
    allSeqs->idArray = realloc(allSeqs->idArray, sizeof(char*) * (allSeqs->numSeqs + 1));
    allSeqs->seqArray = realloc(allSeqs->seqArray, sizeof(char*) * (allSeqs->numSeqs + 1));

    //copy over the ID and the seq
    allSeqs->idArray[allSeqs->numSeqs] = idCopy;
    allSeqs->seqArray[allSeqs->numSeqs] = seqCopy;

    allSeqs->numSeqs++;

}

