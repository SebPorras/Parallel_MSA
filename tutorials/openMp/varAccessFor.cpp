#include <iostream>
#include <cmath>
#include <omp.h> 
#include <chrono> 

int main(int argc, char** argv) 
{
    #include <iostream>
#include <cmath>
#include <omp.h> 
#include <chrono> 


    //how many threads we will use 
    omp_set_num_threads(4); 

    /*
    int i = 10; 

    #pragma omp parallel firstprivate(i) //will use first value of i going forward
    { 
        printf("Thread %d says i = %d\n", i, omp_get_thread_num(), i);
        i = 1000
    }    
      
    printf("actually i = %d\n", i);
    */

   /*
    int i, j; 

    #pragma omp parallel private(i, j) //iterate through their own unique i and j
    { 
        for (i = 0; i < 4; ++i) {
            for (j = 0; j < 10; ++j) {
                printf("Threads %d got i = %d j = %d", omp_get_thread_num(), i, j);
            }
        }
    }  
    */  

  
   /*
    int i, j; 
    #pragma omp parallel private(i, j) { 
        //schedule determines how many chunks are passed to each thread at 1 time
        #pragma omp for  schedule(static, 1)//splits up the work amongst all threads, 1 thread gets each i 
        for (i = 0; i < 4; ++i) { //only outer loop is parallelised 
            for (j = 0; j < 10; ++j) {
                printf("Threads %d got i = %d j = %d", omp_get_thread_num(), i, j);
            }
        }
    }
    */

   /*
    int i, j; 
    #pragma omp parallel private(i, j) { 
        //this will make 2 threads sitting there doing nothing - experiment 
        #pragma omp for  schedule(static, 2)//splits up the work amongst all threads, 1 thread gets each i 
        for (i = 0; i < 4; ++i) { //only outer loop is parallelised 
            for (j = 0; j < 10; ++j) {
                printf("Threads %d got i = %d j = %d", omp_get_thread_num(), i, j);
            }
        }
    }
    */

    /*
    int i, j; 
    #pragma omp parallel private(i, j) { 
        //assign work as it becomes available 
        #pragma omp for  schedule(dynamic)//splits up the work amongst all threads, 1 thread gets each i 
        for (i = 0; i < 4; ++i) { //only outer loop is parallelised 
            for (j = 0; j < 10; ++j) {
                printf("Threads %d got i = %d j = %d", omp_get_thread_num(), i, j);
            }
        }
    }
    */

   /*
    int i, j; 
    #pragma omp parallel private(i, j) { 
        //choose layer number - will iterate over pairs of values 
        #pragma omp for collapse(2) schedule(dynamic) 
        for (i = 0; i < 4; ++i) { 
            for (j = 0; j < 10; ++j) {
                printf("Threads %d got i = %d j = %d", omp_get_thread_num(), i, j);
            }
        }
    }
    */

    int i, j; 

    //once parallel is done, assigns last values of i and j back to original vars
    #pragma omp parallel  { 

        #pragma omp for lastprivate(i, j)
        for (i = 0; i < 4; ++i) { 
            for (j = 0; j < 10; ++j) {
                printf("Threads %d got i = %d j = %d", omp_get_thread_num(), i, j);
            }
        }
    }

    printf("In the end i=%d, j=%d",i , j);

}