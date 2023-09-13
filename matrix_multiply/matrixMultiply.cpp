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
  
    memset(C, 0.0f, N * N * sizeof(float));

    const int ROLL_FACTOR = 4; 

    omp_set_num_threads(32);

    for (int i = 0; i < N; i += 8 * ROLL_FACTOR) {   

        #pragma omp parallel for 
        for (int j = 0; j < N; j++) {

            __m256 c1 = _mm256_load_ps(&C[i + j * N]);
            __m256 c2 = _mm256_load_ps(&C[i + 8 + j * N]);
            __m256 c3 = _mm256_load_ps(&C[i + 16 + j * N]);
            __m256 c4 = _mm256_load_ps(&C[i + + 24 + j * N]);

            
            for (int k = 0; k < N; k++) {

                __m256 b = _mm256_broadcast_ss(&B[j * N + k]);
                
                c1 = _mm256_add_ps(c1,
                 _mm256_mul_ps(_mm256_load_ps(&A[k * N + i]), b));
                
                c2 = _mm256_add_ps(c2,
                 _mm256_mul_ps(_mm256_load_ps(&A[k * N + 8 + i]), b));

                c3 = _mm256_add_ps(c3,
                 _mm256_mul_ps(_mm256_load_ps(&A[k * N + 16 +  i]), b));

                c4 = _mm256_add_ps(c4,
                 _mm256_mul_ps(_mm256_load_ps(&A[k * N + 24 + i]), b));
            }    

            _mm256_store_ps(&C[i + j * N], c1);

            _mm256_store_ps(&C[i + 8 + j * N], c2);

            _mm256_store_ps(&C[i + 16 + j * N], c3);

            _mm256_store_ps(&C[i + 24 + j * N], c4);
        }
    }
}