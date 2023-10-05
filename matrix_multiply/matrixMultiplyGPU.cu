#include "matrixMultiplyGPU.cuh"

__host__ void matrixMultiply_GPU(int N, const float* A, const float* B, float* C, int *arg, int argCount)
{
    int M_LEN = N * N; 

    float* d_A; 
    float* d_B; 
    float* d_C; 

    //allocate memory on the device 
    cudaMalloc(&d_A, sizeof(float) * M_LEN);
    cudaMalloc(&d_B, sizeof(float) * M_LEN); 
    cudaMalloc(&d_C, sizeof(float) * M_LEN);

    //copy across arrays to the device 
    cudaMemcpy(d_A, A, sizeof(float) * M_LEN, cudaMemcpyHostToDevice); 
    cudaMemcpy(d_B, B, sizeof(float) * M_LEN, cudaMemcpyHostToDevice); 
    cudaMemcpy(d_C, C, sizeof(float) * M_LEN, cudaMemcpyHostToDevice); 

    int NUM_THREADS = 32; //32 * 32 is 1024 which is upper limit for threads
    //essentially one block per row and column  
    int NUM_BLOCKS = (N + NUM_THREADS - 1) / NUM_THREADS; 

    dim3 grid(NUM_BLOCKS, NUM_BLOCKS); //use a 2D grid representation 
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

    //make sure we don't go over edge of memory 
    if ((row < N) && (col < N)) {
        float sum = 0.0f; 
        //each thread handles a single row and column multiplication 
        for (int k = 0; k < N; ++k) {
            sum += A[k * N + col] * B[row * N + k]; 
        }

        C[row * N + col] = sum; 
    }
}