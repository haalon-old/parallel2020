#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "problem.h"
#include "block.hpp"

// struct Communicator
// {
//     int rank, size;
//     static Communicator ** array;
//     Communicator(int rank, int size) {
//         this->rank = rank;
//         this->size = size;

//         if(!rank)
//             array = new Communicator*[size];

//         array[rank] = this;
//     }

//     ~Communicator() {
//         if(!rank)
//             delete[] array;
//     }

//     void send(int to, double * buff) {
//         return;
//     }
    
//     void recv(int from, double * buff) {

//     }
// };

void loop()
{
    printf("N=%d K=%d H=%.6f TAU=%.6f\n", N, K, H_X, TAU);
    printf("C_X=%.6f C_Y=%.6f C_Z=%.6f\n", C_X, C_Y, C_Z);

    Block b = Block(0);
    b.init0();
    printf("Error #%3d: %.17f\n", b.t, b.get_error());
    b.init1();
    
    printf("Error #%3d: %.17f\n", b.t, b.get_error());


    for(int t=2; t<=K; t++)
    {
        b.calcNext();
        printf("New Error #%3d: %.17f\n", t, b.get_error());
    }
}

int main(int argc, char const *argv[])
{
    double start = omp_get_wtime();
    loop();
    printf("%f s elapsed\n", (omp_get_wtime() - start));

    return 0;
}
