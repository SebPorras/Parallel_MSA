#include "matrixMultiply.h"
/**
* @brief Implements an NxN matrix multiply C=A*B
*
* @param[in] N : dimension of square matrix (NxN)
* @param[in] A : pointer to input NxN matrix
* @param[in] B : pointer to input NxN matrix
* @param[out] C : pointer to output NxN matrix
* @param[in] args : pointer to array of integers which can be used for debugging and performance tweaks. Optional. If unused, set to zero
* @param[in] argCount : the length of the flags array
* @return void
* */
void matrixMultiply(int N, const float* A, const float* B, float* C, int* args, int argCount)
{
    const int bk = args[0];
    const int bj = args[1]; 
    const int bi = args[2];

    for (int i = 0; i < N; i += bi) {
        for (int j = 0; j < N; j += bj) {
            for (int k = 0; k < N; k += bk) {

                for (int ii = 0; ii < bi; ii++) {
                    
                    int I = ii + i;
                    
                    for (int jj = 0; jj < bj; jj++) {
                        int J = jj + j; 
                        
                        for (int kk = 0; kk < bk; kk++) {
                            int K = k + kk; 
                            C[I * N + J] += A[I * N +K] + B[K * N + J];
                        }

                    }
                }
            }

        }
    } 

}