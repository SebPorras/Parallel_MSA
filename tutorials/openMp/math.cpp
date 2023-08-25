#include <iostream>
#include <cmath>
#include <omp.h> 
#include <chrono> 
//compile command g++ OMPsandbox.cpp -O3 -Wall -Wpedantic -fopenmp 
int main(int argc, char** argv) {

    /*
    //how many threads we will use 
    omp_set_num_threads(4); 
    for (int iter = 0; iter < 16; iter++) {
        
        int count = 0; 

        #pragma omp parallel for 
        for (int i = 1; i <= 40000; i++) {
            
            //very slow 
            #pragma omp critical //only one thread will operate at once
            count += i; 
        }

        printf("Iter %d, count=%d\n", iter, count);
    }
    */

    /*
    //how many threads we will use 
    omp_set_num_threads(4); 
    for (int iter = 0; iter < 16; iter++) {
        
        int count = 0; 

        #pragma omp parallel for 
        for (int i = 1; i <= 40000; i++) {
        
            //a particular update of a var is the only thing that is critical 
            #pragma omp atomic update 
            count += i; 
        }

        printf("Iter %d, count=%d\n", iter, count);
    }
    */

    /*
    omp_set_num_threads(4); 
    for (int iter = 0; iter < 16; iter++) {
        
        int count; 
        int partial_count; 

        //partial is private to each thread but count is shared bewteen 
        #pragma omp parallel private(partial_count) shared(count)
        {
            count = 0;
            partial_count = 0; 

            #pragma omp for 
            for (int i = 1; i <= 40000; i++) {
                //a particular update of a var is the only thing that is critical 
                partial_count += i; 
            }

            #pragma omp critical 
            count += partial_count; 
        }

        printf("Iter %d, count=%d\n", iter, count);
    }
    */

    /*
    omp_set_num_threads(4); 
    for (int iter = 0; iter < 16; iter++) {
        
        int count; 
        int partial_count; 

        //partial is private to each thread but count is shared bewteen 
        #pragma omp parallel private(partial_count) shared(count)
        {
            count = 0;
            partial_count = 0; 

            #pragma omp for 
            for (int i = 1; i <= 40000; i++) {
                //a particular update of a var is the only thing that is critical 
                partial_count += i; 
            }

            #pragma omp critical 
            count += partial_count; 

            //all threads have to update count within an iteration before printing 
            #pragma omp barrier //prevent race conditoin to show count 
            printf("Iter %d, count=%d\n", iter, count);
        }       
    }
    */

   /*
    omp_set_num_threads(4); 
    for (int iter = 0; iter < 16; iter++) {
        
        int count; 
        int partial_count; 

        //partial is private to each thread but count is shared bewteen 
        #pragma omp parallel private(partial_count) shared(count)
        {
            count = 0;
            partial_count = 0; 

            #pragma omp for 
            for (int i = 1; i <= 40000; i++) {
                //a particular update of a var is the only thing that is critical 
                partial_count += i; 
            }

            #pragma omp critical 
            count += partial_count; 

            //all threads have to update count within an iteration before printing 
            #pragma omp barrier //prevent race conditoin to show count 
            #pragma omp single 
            printf("Iter %d, thread=%d count=%d\n", iter, get_comp_thread_num(), count);
        }       
    }
    */

   /*
   omp_set_num_threads(4); 
   for (int iter = 0; iter < 16; iter++) {
        
        int count; 
        int partial_count; 

        //partial is private to each thread but count is shared bewteen 
        #pragma omp parallel private(partial_count) shared(count)
        {
            count = 0;
            partial_count = 0; 

            #pragma omp for 
            for (int i = 1; i <= 40000; i++) { 
                partial_count += i; 
            }

            #pragma omp critical 
            count += partial_count; 
 
            #pragma omp barrier 
            //the base thread will handle the work 
            #pragma omp master
            printf("Iter %d, thread=%d count=%d\n", iter, get_comp_thread_num(), count);
        }       
    }
    */

    //downside is that threads will sit idle if you assign more threads then sections
    omp_set_num_threads(4); 
    for (int iter = 0; iter < 5; iter++) {

        //can split up segments of code that aren't just for loops 
        #pragma omp parallel sections 
        { 
            #pragma omp section 
            {
                printf("Thread %d is in section 1\n", omp_get_thread_num());
                //did some math 
            }

            #pragma omp section 
            {
                printf("Thread %d is in section 2\n", omp_get_thread_num());
                //did some other math 
            }

            printf("--------------------------------------\n");

        }


    }
  

}