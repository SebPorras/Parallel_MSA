#include <chrono>
#include <cstdio>
#include <cmath>
#include <mpi.h> 

// Simple sum over an array.
double sum(const double* f, int N) {
    double s = 0.;
    for (int i = 0; i < N; ++i)
        s += f[i];
    return s;
}

int main(int argc, char **argv) {

    //initialise the MPI program 
    MPI_Init(&argc, &argv);

    MPI_Status status;

    //work out how many processes are running 
    //a communicator is what manages the task 
    int worldSize; 
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    //for each task what number are you 
    int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("I'm process %d of %d\n", rank, worldSize);  

    const double lower = 0.;
    const double upper = 1.;
    const int N = 1000000;
    const double dx = (upper - lower) / N;

    //define some unique pieces based on rank and world size 
    int NPerRank = int(float(N) / float(worldSize)); 
    int myFirstN = rank * NPerRank;
    int myLastN = (rank + 1) * NPerRank; 
    //last process is assigned the last chunk of work 
    if (my_rank == worldSize - 1) myLastN = N; 

    //make sure we compute a unique part of the integral. 
    double myLower = myFirstN * dx; 

    // Populate f(xi) can split N across multiple processes 
    //can lower memory usage because we do a fraction of N 
    //memory is not shared by MPI processes so you save on memory usage 
    double* f = (double*)malloc(sizeof(*f) * NPerRank);
    for (int i = 0; i < NPerRank; ++i) {
        double xi = myLower + (i + 0.5) * dx;
        f[i] = 1. / (xi * xi + 1.);
    }

    // Begin timing
    auto start = std::chrono::high_resolution_clock::now();

    // Calculate integral
    double integral = sum(f, NPerRank) * dx;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> duration_us = end - start;

    double result = integral; 

    /*
    if (rank == 0) {
        //get ready to recieve all the pieces of the integral from 
        // other processes. 
        double other_piece; //will store portion of other results 
        for (int neighbour = 1; neighbour < worldSize; ++neighbour) {
            MPI_Recv(&other_piece, 1, MPI_DOUBLE, neighbour, 
            100 + neighbour, MPI_COMM_WORLD, &status); 
            result += other_piece; 
        }

        printf("Integral was %.15g\n", result);
        printf("Error was %.10e\n", std::abs(result - M_PI/4.));
        printf("Time taken: %g us\n", duration_us.count());
        
    } else {
        //pass my result to the master process. 
        //data to send, things to send, data size, which rank to send to,
        //identifier tag, communicator)
        MPI_Send(&integral, 1, MPI_DOUBLE, 0, 100 + rank, MPI_COMM_WORLD); 
    }
    Can reduce this code above to a single line using MPI reduce 
    */

    /*
    //(thing we're reducing on, where to store result, how many things are being recieved,
    //data type, the operation to perform, which process to return to, comm world)
    MPI_Reduce(&integral, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 

    //also don't have to use the broadcast if we use MPI_allreduce,
    // combines reduce and broadcast
    //need to broadcast integral back to other processes
    MPI_Bcast(&result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);    
    */ 
    MPI_Allreduce(&integral, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    //now do other things with the value of teh integral 
    printf("Process %d thinks integral is %lf\n", rank, result); 
 


    free(f);

    MPI_Finalize(); 
}
