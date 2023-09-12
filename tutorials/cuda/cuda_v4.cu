#include <chrono>
#include <cstdio>
#include <cmath>

// This macro lets us output the file and line number
// if a CUDA error occurs.
#define cudaCheck(expr) \
    do { \
        cudaError_t e = (expr); \
        if (e != cudaSuccess) { \
            fprintf(stderr, "CUDA error: %s (%s:%d)\n", cudaGetErrorString(e), __FILE__, __LINE__); \
            abort(); \
        } \
    } while (false)


// Works with any number of blocks.
// Requires shared memory space for num_threads doubles
// Requires f to be aligned as double4 (256 bit)
__global__
void sum(const double* f, int N, double* out) {
    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Serial sum using vectorised memory loads
    double s = 0.;
    int i = idx;
    for (; i < N/4; i += gridDim.x * blockDim.x) {
        // Safe since f is 256bit aligned (sizeof(double4) == 256bit)
        double4 tmp = reinterpret_cast<const double4*>(f)[i];
        s += tmp.x + tmp.y + tmp.z + tmp.w;
    }

    // Handle left over elements if N % 4 != 0
    if (4*i < N) {
        double4 tmp = reinterpret_cast<const double4*>(f)[i];
        switch (N - 4*i) {
            case 3:
                s += tmp.z;
            case 2:
                s += tmp.y;
            case 1:
                s += tmp.x;
        }
    }

    // Store result of serial sum in shared memory
    extern __shared__ double smem[];
    smem[tid] = s;

    // Perform reduction over shared memory
    for (int i = blockDim.x/2; i > 0; i /= 2) {
        __syncthreads();
        if (tid < i)
            smem[tid] += smem[tid + i];
    }
    if (tid == 0) out[blockIdx.x] = smem[0];
}

int main() {
    const double lower = 0.;
    const double upper = 1.;
    const int N = 1'000'000;

    const double dx = (upper - lower) / N;

    // Populate f(xi)
    double* f = (double*)malloc(sizeof(*f) * N);
    for (int i = 0; i < N; ++i) {
        double xi = lower + (i + 0.5) * dx;
        f[i] = 1. / (xi * xi + 1.);
    }

    // Enable 8 byte shared memory banks for conflict-free access to consecutive doubles
    cudaCheck(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));

    // Set up working memory
    const int NUM_THREADS = 1024;
    int num_blocks = ((N+3)/4 + NUM_THREADS - 1) / NUM_THREADS; // Only need N/4 threads

    double* d_f; cudaCheck(cudaMalloc((void**)&d_f, sizeof(*d_f) * N));
    // Over-allocating the result array to leave space for multiple reduction steps
    double* d_result; cudaCheck(cudaMalloc((void**)&d_result, sizeof(*d_result) * num_blocks * 2));

    // Copy f(xi) array to GPU
    cudaCheck(cudaMemcpy(d_f, f, sizeof(*f) * N, cudaMemcpyHostToDevice));

    // Begin timing
    auto start = std::chrono::high_resolution_clock::now();

    // Calculate integral
    sum<<<num_blocks, NUM_THREADS, sizeof(*f)*NUM_THREADS>>>(d_f, N, d_result);
    double* buf1 = d_result;
    double* buf2 = &d_result[num_blocks];
    while (num_blocks > 1) {
        int old_blocks = num_blocks; // This will be the number of elements we need to reduce over
        num_blocks = ((num_blocks+3)/4 + NUM_THREADS - 1) / NUM_THREADS;
        sum<<<num_blocks, NUM_THREADS, sizeof(*f)*NUM_THREADS>>>(buf1, old_blocks, buf2);
        double* tmp = buf2;
        buf2 = buf1;
        buf1 = tmp;
    }

    // Copy final result back to the host
    double integral;
    cudaCheck(cudaMemcpy(&integral, buf1, sizeof(integral), cudaMemcpyDeviceToHost));
    integral *= dx;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration_us = end - start;

    printf("Integral was %.15g\n", integral);
    printf("Error was %.10e\n", std::abs(integral - M_PI/4.));
    printf("Time taken: %g us\n", duration_us.count());

    cudaFree(d_f);
    cudaFree(d_result);
    free(f);
}
