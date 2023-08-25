#include <iostream>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <random>

int main(int argc, char **argv) 
{

	// Let's evaluate e^PI using Monte Carlo techniques and openMP
    printf("|| Nsamples  |   Value   |   Error  | Time (us)||\n");

    for (int niter=3; niter<9; niter++) {

    	// Begin timing
		auto start = std::chrono::high_resolution_clock::now();

	    int Nsamples = int(pow(10, niter));

		// We can evaluate 'PI' by drawing two random numbers x and y from 0 to 1 (i.e., throwing darts at a uniform square)
		// We then count how many pairs of x and y satisfy x^2 + y^2 = 1 (i.e., fall within a quarter of a circle with radius 1 inside the square)
		// The ratio of the number of pairs inside the circle vs. the total number of pairs is equal to the area of the quarter circle (PI/4) divided
		// by the area of the square (1). And so this ratio must be PI/4.

		int count = 0;
        int total = 0;

        #pragma omp parallel  
        {
            // First, set up a random number generator
            thread_local std::mt19937 gen(std::random_device{}());
            //every thread needs its own version of RNG (thread safe version)
            thread_local std::uniform_real_distribution<> dist(0, 1);

            #pragma omp for reduction(+:count)
            for (int i = 0; i < Nsamples; i++) {

                // Generate two random numbers between 0 and 1
                double x = dist(gen);
                double y = dist(gen);

                // Compute radius squared and add to the counter if it's within the unit circle.
                double r2 = x*x + y*y;
                if (r2 <= 1.0) count++;

            }

            // We can evaluate 'e' almost by magic by drawing random numbers from 0 to 1 and then adding them sequentially until the total is larger than 1.
            // Once it is larger than 1, we add how many values we needed to draw to reach that point to a list of samples, then take their mean. 
            // The mean is e (Trust me, it works, because fundamentally e = Sum(1/n!). But it seems like magic...)

            #pragma omp for reduction(+:total)
            for (int i = 0; i < Nsamples; i++) {

                int count = 0;
                double sample = 0.0;

                // Add a random number to the total between 0 and 1 until it is larger than 1.
                while(sample <= 1.0) {
                    count++;
                    sample += dist(gen);
                }
                
                // Add it to the running total
                total += count;
            }
        }
		
	    double ratio = double(count)/double(Nsamples);
	    double eval = double(total)/double(Nsamples);
	    double answer = pow(eval,4.0*ratio);

	    auto end = std::chrono::high_resolution_clock::now();
	    std::chrono::duration<double, std::micro> duration_us = end - start;

	    printf("|| %9d | %3.6lf | %3.6lf | %8.0lf ||\n", Nsamples, answer, std::abs(answer - exp(M_PI)), duration_us.count());
	}

}