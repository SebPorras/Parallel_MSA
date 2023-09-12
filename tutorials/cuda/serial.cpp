#include <chrono>
#include <cstdio>
#include <cmath>
#include <stdlib.h>

// Simple sum over an array.
double sum(const double* f, int N) {
    double s = 0.;
    for (int i = 0; i < N; ++i)
        s += f[i];
    return s;
}

int main() {
    const double lower = 0.;
    const double upper = 1.;
    const int N = 1000000;

    const double dx = (upper - lower) / N;

    // Populate f(xi)
    double* f = (double*)malloc(sizeof(*f) * N);
    for (int i = 0; i < N; ++i) {
        double xi = lower + (i + 0.5) * dx;
        f[i] = 1. / (xi * xi + 1.);
    }

    // Begin timing
    auto start = std::chrono::high_resolution_clock::now();

    // Calculate integral
    double integral = sum(f, N) * dx;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration_us = end - start;

    printf("Integral was %.15g\n", integral);
    printf("Error was %.10e\n", std::abs(integral - M_PI/4.));
    printf("Time taken: %g us\n", duration_us.count());

    free(f);
}
