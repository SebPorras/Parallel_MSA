#include "matrixMultiplyMPI.h"

void matrixMultiply_MPI(int N, const float* A, const float* B, 
                        float* C, int* args, int argCount) {

    int worldSize; 
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int NPerRank = int(float(N * N) / float(worldSize)); //(2048 * 2048) / 2

    int firstN = (rank  == 0) ? 0 : NPerRank; 
    int lastN = (rank == 0) ? NPerRank : N * N; //2048 * 2048

    float* localC = (float*) malloc(sizeof(float) * (N * N)); 
    float* localA = (float*) malloc(sizeof(float) * NPerRank);
    for (int i = 0; i < NPerRank; ++i) {
        localA[i] = A[i + firstN]; 
    }

    matrixMultiply(NPerRank, localA, B, localC, args, argCount); 
    MPI_Allgather(localC, NPerRank, MPI_FLOAT, C, NPerRank,
     MPI_FLOAT, MPI_COMM_WORLD); 


    free(localA); 
    free(localC); 
}

