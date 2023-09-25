#include "matrixMultiplyGPU.cuh"
#include <stdlib.h>
#include <stdio.h>

#define cudaCheck(expr) \
    do { \
        cudaError_t e = (expr); \
        if (e != cudaSuccess) { \
            fprintf(stderr, "CUDA error: %s (%s:%d)\n", cudaGetErrorString(e), __FILE__, __LINE__); \
            abort(); \
        } \
    } while (false)



__host__ void matrixMultiply_GPU(int N, const float* A, const float* B, float* C, int *arg, int argCount)
{
    memset(C, 0.0f, N * N * sizeof(float)); 
    
    int M_LEN = N * N; 

    float* d_A; 
    float* d_B; 
    float* d_C; 

    cudaCheck(cudaMalloc(&d_A, sizeof(float) * M_LEN));
    cudaCheck(cudaMalloc(&d_B, sizeof(float) * M_LEN)); 
    cudaCheck(cudaMalloc(&d_C, sizeof(float) * M_LEN)); 

    cudaCheck(cudaMemcpy(d_A, A, sizeof(float) * M_LEN, cudaMemcpyHostToDevice)); 
    cudaCheck(cudaMemcpy(d_B, B, sizeof(float) * M_LEN, cudaMemcpyHostToDevice)); 
    cudaCheck(cudaMemcpy(d_C, C, sizeof(float) * M_LEN, cudaMemcpyHostToDevice)); 

    int NUM_THREADS = 32; 
    int NUM_BLOCKS = N / NUM_THREADS; 

    dim3 grid(NUM_BLOCKS, NUM_BLOCKS); 
    dim3 threads(NUM_THREADS, NUM_THREADS); 
    
    matrixMultiplyKernel_GPU<<<grid, threads>>>(N, d_A, d_B, d_C, 0, 0, 0); 

    cudaCheck(cudaMemcpy(C, d_C, sizeof(float) * M_LEN, cudaMemcpyDeviceToHost)); 

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C); 	
}

__global__ void matrixMultiplyKernel_GPU(int N, const float* A, const float* B, float* C, int flag0, int flag1, int flag2)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;  

    float sum = 0.0f; 

    if ((row < N) && (col < N)) {
        
        for (int k = 0; k < N; ++k) {
            sum += A[row * N + k] * B[k * N + col]; 
        }

        C[row * N + col] = sum; 
    }
}

