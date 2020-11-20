#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "problem.h"
#include "block.hpp"
#include "comm.hpp"


void loop()
{
    printf("N=%d K=%d H=%.6f TAU=%.6f\n", N, K, H_X, TAU);
    printf("C_X=%.6f C_Y=%.6f C_Z=%.6f\n", C_X, C_Y, C_Z);

    Block b = Block(0);
    Comm c = Comm(0, 1, &b);

    b.init0();
    printf("Error #%3d: %.17f\n", b.t, b.get_error());
    b.init1(&c);
    
    printf("Error #%3d: %.17f\n", b.t, b.get_error());


    for(int t=2; t<=K; t++)
    {
        b.calcNext(&c);
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
