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
// `out` should have space for num_blocks values
// `buf` should have space for num_blocks * num_threads values
__global__
void sum(const double* f, int N, double* out, double* buf) {
    int tid = threadIdx.x;
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // serial sum with a stride of grid size
    double s = 0.;
    for (int i = idx; i < N; i += gridDim.x * blockDim.x)
        s += f[i];

    // Store result of serial sum in buffer
    buf[idx] = s;

    // Perform reduction over this block, using the buffer
    for (int i = blockDim.x / 2; i > 0; i /= 2) {
        __syncthreads();
        if (tid < i)
            buf[idx] += buf[idx + i];
    }

    // Store result for each block in output
    if (tid == 0) out[blockIdx.x] = buf[idx];
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


    // Set up working memory
    const int NUM_THREADS = 1024;
    int num_blocks = (N + NUM_THREADS - 1) / NUM_THREADS;
    double* d_f; cudaCheck(cudaMalloc((void**)&d_f, sizeof(*d_f) * N));
    double* d_buffer; cudaCheck(cudaMalloc((void**)&d_buffer, sizeof(*d_buffer) * num_blocks * NUM_THREADS));
    // Over-allocating the result array to leave space for multiple reduction steps
    double* d_result; cudaCheck(cudaMalloc((void**)&d_result, sizeof(*d_result) * num_blocks * 2));

    // Copy f(xi) array to GPU
    cudaCheck(cudaMemcpy(d_f, f, sizeof(*f) * N, cudaMemcpyHostToDevice));

    // Begin timing
    auto start = std::chrono::high_resolution_clock::now();

    // Calculate integral

    // Initial reduction will give num_blocks values in d_result
    sum<<<num_blocks, NUM_THREADS>>>(d_f, N, d_result, d_buffer);

    // Keep applying reduction to d_result until only one block is needed
    // (so the result of that block will be the full sum of f).
    // Uses extra space allocated in d_result to store the result of
    // the next reduction step (buf2).
    double* buf1 = d_result;
    double* buf2 = &d_result[((num_blocks+3)/4) * 4]; // Make sure buf2 is aligned to 256 bits
    while (num_blocks > 1) {
        int old_blocks = num_blocks; // This will be the number of elements we need to reduce over
        num_blocks = (num_blocks + NUM_THREADS - 1) / NUM_THREADS;
        sum<<<num_blocks, NUM_THREADS>>>(buf1, old_blocks, buf2, d_buffer);

        // The result is in buf2. We want to reduce over that result, so swap
        // the pointers to make buf1 point to it, and buf2 now points to usable
        // memory for the next result.
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
    cudaFree(d_buffer);
    free(f);
}
