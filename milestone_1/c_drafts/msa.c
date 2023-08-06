//#include "hdf5.h"
#include "io.h"

int main(int argc, char** argv) {

    //check the user provdes a file 
    if (argc == 1) {
        fprintf(stderr, "Provide a fasta file\n");
        return 1;
    }

    FILE* fileStream;
    fileStream = fopen(argv[FILENAME], "r");

    if (!fileStream) {
        fprintf(stderr, "File read error\n");
        return 2;
    }

    FastaSeqs seqs = read_fasta_file(fileStream);
    fclose(fileStream);

    for (int i = 0; i < seqs.numSeqs; ++i) {
        printf("%s %s\n", seqs.idArray[i], seqs.seqArray[i]);
    }

    printf("%d\n", seqs.numSeqs);

    free(seqs.idArray);
    free(seqs.seqArray);
    
    return 0;
}