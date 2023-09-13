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

    const int ROLL_FACTOR = 8; 
    const int VEC_LEN = 8; 

    omp_set_num_threads(64);

    for (int i = 0; i < N; i += VEC_LEN * ROLL_FACTOR) {   

        #pragma omp parallel for 
        for (int j = 0; j < N; j++) {

            __m256 output_vecs[ROLL_FACTOR]; 

            for (int x = 0; x < ROLL_FACTOR; x += 4) {
                output_vecs[x] = _mm256_load_ps(&C[i + (VEC_LEN * x) + j * N]);
                output_vecs[x + 1] = _mm256_load_ps(&C[i + (VEC_LEN * (x + 1)) + j * N]);
                output_vecs[x + 2] = _mm256_load_ps(&C[i + (VEC_LEN * (x + 2)) + j * N]);
                output_vecs[x + 3] = _mm256_load_ps(&C[i + (VEC_LEN * (x + 3)) + j * N]);
            }
                  
            for (int k = 0; k < N; k++) {

                __m256 b = _mm256_broadcast_ss(&B[j * N + k]);
                
                for (int x = 0; x < ROLL_FACTOR; x += 4) {
                    output_vecs[x] = _mm256_add_ps(output_vecs[x],
                                     _mm256_mul_ps(
                                     _mm256_load_ps(
                                        &A[k * N + (VEC_LEN * x) + i]), b));

                    output_vecs[x + 1] = _mm256_add_ps(output_vecs[x + 1],
                                     _mm256_mul_ps(
                                     _mm256_load_ps(
                                        &A[k * N + (VEC_LEN * (x + 1)) + i]), b));

                    output_vecs[x + 2] = _mm256_add_ps(output_vecs[x + 2],
                                     _mm256_mul_ps(
                                     _mm256_load_ps(
                                        &A[k * N + (VEC_LEN * (x + 2)) + i]), b));

                    output_vecs[x + 3] = _mm256_add_ps(output_vecs[x + 3],
                                     _mm256_mul_ps(
                                     _mm256_load_ps(
                                        &A[k * N + (VEC_LEN * (x + 3)) + i]), b));
                }   
            }    

            for (int x = 0; x < ROLL_FACTOR; x += 4) {
                _mm256_store_ps(&C[i + (x * VEC_LEN) + j * N], output_vecs[x]);
                _mm256_store_ps(&C[i + ((x + 1) * VEC_LEN) + j * N], output_vecs[x + 1]);
                _mm256_store_ps(&C[i + ((x + 2) * VEC_LEN) + j * N], output_vecs[x + 2]);
                _mm256_store_ps(&C[i + ((x + 3) * VEC_LEN) + j * N], output_vecs[x + 3]);
            }
        }
    }
}