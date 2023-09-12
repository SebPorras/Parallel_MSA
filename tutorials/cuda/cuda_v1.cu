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


// Only works with a single block!
__global__
void sum(const double* f, int N, double* out) {
    int idx = threadIdx.x;

    // serial sum with a stride of the block size
    double s = 0.;
    for (int i = idx; i < N; i += blockDim.x)
        s += f[i];

    // Store result of serial sum in buffer
    out[idx] = s;

    // Perform reduction over buffer
    for (int i = blockDim.x / 2; i > 0; i /= 2) {
        __syncthreads(); // Need to synchronize all threads in the block between steps
        if (idx < i)
            out[idx] += out[idx + i];
    }

    // Now out[0] contains the full sum.
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
    double* d_result; cudaCheck(cudaMalloc((void**)&d_result, sizeof(*d_result) * NUM_THREADS));
    double* d_f; cudaCheck(cudaMalloc((void**)&d_f, sizeof(*d_f) * N));

    // Copy f(xi) array to GPU
    cudaCheck( (d_f, f, sizeof(*f) * N, cudaMemcpyHostToDevice));

    // Begin timing
    auto start = std::chrono::high_resolution_clock::now();

    // Calculate integral
    sum<<<1, NUM_THREADS>>>(d_f, N, d_result);

    // Copy result back to host
    double integral;
    cudaCheck(cudaMemcpy(&integral, d_result, sizeof(integral), cudaMemcpyDeviceToHost));
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
