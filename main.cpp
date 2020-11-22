#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "problem.h"
#include "block.hpp"
#include "comm.hpp"

void sync(Block ** b, Comm ** c)  {
    float m = 0;
    // for(int i = 0; i<BX*BY*BZ; i++)
    //     b[i]->prepare();
    for(int i = 0; i<BX*BY*BZ; i++)
    {
        b[i]->exchange(c[i]);
        m = b[i]->get_error() > m ? b[i]->get_error() : m;
    }
    printf("\tGlobal err: %.17f\n", m);
}

void loop()
{
    printf("N=%d K=%d H=%.6f TAU=%.6f\n", N, K, H_X, TAU);
    printf("C_X=%.6f C_Y=%.6f C_Z=%.6f\n", C_X, C_Y, C_Z);
    Block * b[BX*BY*BZ];
    Comm * c[BX*BY*BZ];

    printf("Error #%3d\n", 0);
    for(int i = 0; i<BX*BY*BZ; i++) {
        b[i] = new Block(i);
        c[i] = new Comm(i, BX*BY*BZ, b[i]);

        b[i]->init0();
        printf("\tBlock #%3d: %.17f\n", i, b[i]->get_error());
    }
    sync(b,c);
    printf("\n");

    printf("Error #%3d\n", 1);
    for(int i = 0; i<BX*BY*BZ; i++) {
        b[i]->init1();
        printf("\tBlock #%3d: %.17f\n", i, b[i]->get_error());

    }
    sync(b,c);
    printf("\n");

    for(int t=2; t<=K; t++)
    {
        printf("Error #%3d\n", t);
        for(int i = 0; i<BX*BY*BZ; i++) {
            b[i]->calcNext();
            printf("\tBlock #%3d: %.17f\n", i, b[i]->get_error());
        }
        sync(b,c);
        printf("\n");
    }
}

int main(int argc, char const *argv[])
{
    double start = omp_get_wtime();
    loop();
    printf("%f s elapsed\n", (omp_get_wtime() - start));

    return 0;
}
