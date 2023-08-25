#include <chrono>
#include <cstdio>
#include <cmath>
#include <omp.h>

int main(int argc, char **argv) {

    printf("Running with %d threads\n", omp_get_max_threads());
    // Begin timing
    auto start = std::chrono::high_resolution_clock::now();

    int const num_steps = 100000000;
    double const step = 1.0 / num_steps;
    double sum = 0;

    /* Standard implementation 
    #pragma omp parallel 
    { 
        double local_sum = 0;

        #pragma omp for 
        for (int i = 0; i < num_steps; ++i)
        {
            double xi = (i+0.5)*step;
            double f = 1.0 / (xi*xi+1.0);
            local_sum += f;
        }
        
        #pragma omp critical 
        sum += local_sum;
    }
    */

    /*Reduction method 
    //add to the sum based on reduction 
    #pragma omp parallel for reduction( +:sum) 
    for (int i = 0; i < num_steps; ++i)
    {
        double xi = (i+0.5)*step;
        double f = 1.0 / (xi*xi+1.0);
        sum += f;
    }
    */
    
    //ask the compiler to use SIMD 
    #pragma omp parallel for simd reduction( +:sum) 
    for (int i = 0; i < num_steps; ++i)
    {
        double xi = (i+0.5)*step;
        double f = 1.0 / (xi*xi+1.0);
        sum += f;
    }
     

    double result = step * sum;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration_us = end - start;

    printf("Integral was %.15g\n", result);
    printf("Error was %.10e\n", std::abs(result - M_PI/4.));
    printf("Time taken: %g us\n", duration_us.count());

}