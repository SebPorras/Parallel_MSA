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
void matrixMultiply(int N, const float* A, const float* B, float* C,
 int* args, int argCount)  {
    const int bk = 8;
    const int bj = 8; 
    const int bi = 8;

    memset(C, 0.0f, N * N * sizeof(float));

    for (int i = 0; i < N; i += 8) {   
        for (int j = 0; j < N; j++) {

            __m256 c0 = _mm256_setzero_ps();
            for (int k = 0; k < N; k++) {
                
                c0 = _mm256_add_ps(c0, _mm256_mul_ps(
                    _mm256_load_ps(&A[k * N + i]),
                    _mm256_broadcast_ss(&B[j * N + k])));
            }    

            _mm256_store_ps(C + i + j * N, c0); 
        }
    }
}