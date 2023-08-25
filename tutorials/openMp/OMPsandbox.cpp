#include <iostream>
#include <cmath>
#include <omp.h> 
#include <chrono> 

int main(int argc, char** argv) 
{


    //how many threads we will use 
    omp_set_num_threads(8); 

    omp_set_dynamic(1); //tell code at runtime it can use how many threads it thinks is appropriate
    #pragma omp parallel
    {
        for (int i = 0; i < omp_get_num_threads(); i++) 
        {
            if (if == omp_get_thread_num() == i) printf("Thread %d\n", omp_get_num_threads());
            #pragma omp barrier //all threads must reach here before going to next section
            #//will now print in order 
        }
        
    }

}