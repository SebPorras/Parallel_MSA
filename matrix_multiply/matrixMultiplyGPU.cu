#include "matrixMultiplyGPU.cuh"

__host__ void matrixMultiply_GPU(int N, const float* A, const float* B, float* C, int *arg, int argCount)
{
    int M_LEN = N * N; 

    float* d_A; 
    float* d_B; 
    float* d_C; 

    cudaMalloc(&d_A, sizeof(float) * M_LEN);
    cudaMalloc(&d_B, sizeof(float) * M_LEN); 
    cudaMalloc(&d_C, sizeof(float) * M_LEN);
    cudaMemcpy(d_A, A, sizeof(float) * M_LEN, cudaMemcpyHostToDevice); 
    cudaMemcpy(d_B, B, sizeof(float) * M_LEN, cudaMemcpyHostToDevice); 
    cudaMemcpy(d_C, C, sizeof(float) * M_LEN, cudaMemcpyHostToDevice); 

    int NUM_THREADS = 32; 
    int NUM_BLOCKS = (N + NUM_THREADS - 1) / NUM_THREADS; 

    dim3 grid(NUM_BLOCKS, NUM_BLOCKS); 
    dim3 threads(NUM_THREADS, NUM_THREADS); 
    
    matrixMultiplyKernel_GPU<<<grid, threads>>>(N, d_A, d_B, d_C, 0, 0, 0); 

    cudaMemcpy(C, d_C, sizeof(float) * M_LEN, cudaMemcpyDeviceToHost); 

    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C); 	
}

__global__ void matrixMultiplyKernel_GPU(int N, const float* A, const float* B, float* C, int flag0, int flag1, int flag2)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;  

    if ((row < N) && (col < N)) {
        float sum = 0.0f; 
        for (int k = 0; k < N; ++k) {
            sum += A[k * N + col] * B[row * N + k]; 
        }

        C[row * N + col] = sum; 
    }
}