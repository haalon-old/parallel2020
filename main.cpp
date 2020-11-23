#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "problem.h"
#include "block.hpp"
#include "comm.hpp"
#include <string>

#ifndef FAKEMPI
    #include <mpi.h>
#endif

using std::string;

double L_X = 1.0;
double L_Y = 1.0;
double L_Z = 1.0;
double T = 1.0;

int N = 32;
int K = 64;

int BX = 1;
int BY = 1;
int BZ = 1;

double H_X = L_X / N;
double H_Y = L_Y / N;
double H_Z = L_Z / N;

double TAU = T / K;

double C_X = TAU / H_X;
double C_Y = TAU / H_Y;
double C_Z = TAU / H_Z;

// why not
#define SYNCTAG 42

void parseArgs(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i)
    {
        string param(argv[i]);

        if(param == "-n")
            N = atoi(argv[++i]);

        if(param == "-k")
            K = atoi(argv[++i]);

        if(param == "-b") {
            BX = atoi(argv[++i]);
            BY = atoi(argv[++i]);
            BZ = atoi(argv[++i]);
        }

        if(param == "-l") {
            L_X = atof(argv[++i]);
            L_Y = atof(argv[++i]);
            L_Z = atof(argv[++i]);
        }

        if(param == "-t")
            T = atof(argv[++i]);
    }
    H_X = L_X / N;
    H_Y = L_Y / N;
    H_Z = L_Z / N;

    TAU = T / K;

    C_X = TAU / H_X;
    C_Y = TAU / H_Y;
    C_Z = TAU / H_Z;
}
#ifdef FAKEMPI
    void sync(Block ** b, Comm ** c)  {
        float m = 0;
        // for(int i = 0; i<BX*BY*BZ; i++)
        //     b[i]->prepare();
        for(int i = 0; i<BX*BY*BZ; i++)
        {
            c[i]->exchange();
            m = b[i]->get_error() > m ? b[i]->get_error() : m;
        }
        printf("\tGlobal err: %.17f\n", m);
    }

    void loop(int argc, char *argv[]) {
        parseArgs(argc, argv); 
        double start = omp_get_wtime();

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
        printf("%f s elapsed\n", (omp_get_wtime() - start));
    }
#else
    void sync(Block * b, Comm * c) {
        double err = b->get_error(), temp;
        printf("\tBlock #%3d: %.17f\n", c->rank, err);

        c->exchange();

        if(c->rank)
            c->send(0, 1, &err, SYNCTAG);
        else {
            for (int i = 1; i < BX*BY*BZ; ++i) {
                c->recv(i, 1, &temp, SYNCTAG);
                err = temp > err ? temp : err;
            }
            printf("\tGlobal err: %.17f\n", err);
        }
    }

    void loop(int argc, char *argv[]) {
        MPI_Init(&argc, &argv);
        parseArgs(argc, argv); 
        double start = MPI_Wtime();
        

            // Get the number of processes
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // Get the rank of the process
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        Block b(rank);
        Comm c(rank, size, &b);
        if(!rank)
            printf("Error #%3d\n", 0);
        b.init0();
        sync(&b, &c);

        if(!rank)
            printf("Error #%3d\n", 0);
        b.init1();
        sync(&b, &c);

        for(int t=2; t<=K; t++) {
            if(!rank)
                printf("Error #%3d\n", t);
            b.calcNext();
            sync(&b, &c);
        }
        if(!rank)
            printf("%f s elapsed\n", (MPI_Wtime() - start));
        MPI_Finalize();

    }
#endif

int main(int argc, char *argv[]) {   
    loop(argc, argv);    

    return 0;
}
