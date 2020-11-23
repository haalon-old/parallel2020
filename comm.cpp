#include "comm.hpp"
#include "block.hpp"
#include <string.h>
#include <stdio.h>
#include "problem.h"
#ifndef FAKEMPI
    #include <mpi.h>
#endif

int mod(int i, int n) {
    return ((i % n) + n) % n;
}

Comm::Comm() {}

Comm::Comm(int rank, int size, Block * block) {
    this->rank = rank;
    this->block = block;

    int k = rank % BZ;
    int j = rank / BZ % BY;
    int i = rank / BZ / BY;

    neigh[3] = mod(k+1, BZ) + BZ*j + BZ*BY*i;
    neigh[2] = mod(k-1, BZ) + BZ*j + BZ*BY*i;

    neigh[4] = k + BZ*mod(j+1, BY) + BZ*BY*i;
    neigh[1] = k + BZ*mod(j-1, BY) + BZ*BY*i;

    neigh[5] = k + BZ*j + BZ*BY*mod(i+1, BX);
    neigh[0] = k + BZ*j + BZ*BY*mod(i-1, BX);

    // neigh = {mx, my, mz, pz, py, px};
    // printf("\tx %d %d, y %d %d, z %d %d\n", mx, px, my, py, mz, pz);
    #ifdef FAKEMPI
        if(!rank)
            array = new Comm*[size];

        array[rank] = this;
    #endif
}

#ifndef FAKEMPI
    #include <mpi.h>
    Comm::~Comm() {}

    void Comm::send(int to, int size, double * buff, int tag) {
        MPI_Ssend(buff, size, MPI_DOUBLE, to, tag, MPI_COMM_WORLD);
    }

    void Comm::recv(int from, int size, double * buff, int tag) {
        MPI_Recv(buff, size, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    void Comm::exchange() {
        int sizes[3] = {block->nz*block->ny, block->nz*block->nx, block->nx*block->ny};
        for (int i = 0; i <= 2; ++i) {
            int size = sizes[i];
            double * temp1 = new double[size];
            double * temp2 = new double[size];
            double * buff1 = block->edges[i];
            double * buff2 = block->edges[5-i];
            MPI_Request req1;
            MPI_Request req2;
            // MPI_Status status;

            memcpy(temp1, buff1, size*sizeof(double));
            memcpy(temp2, buff2, size*sizeof(double));

            // if temp1 contains edges[0]
            // we want to send it to neigh[0]
            // mark it with tag 5 (coz we want itss edges[5])
            MPI_Isend(temp1, size, MPI_DOUBLE, neigh[i], 5-i, MPI_COMM_WORLD, &req1);
            MPI_Isend(temp2, size, MPI_DOUBLE, neigh[5-i], i, MPI_COMM_WORLD, &req2);

            // we want to recieve edge[0]
            // from neigh[0]
            // it was sent with tag 0
            MPI_Recv(buff1, size, MPI_DOUBLE, neigh[i], i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(buff2, size, MPI_DOUBLE, neigh[5-i], 5-i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Wait(&req1, MPI_STATUS_IGNORE);
            MPI_Wait(&req2, MPI_STATUS_IGNORE);
            delete[] temp1;
            delete[] temp2;
        }
    }
#else
    // define static member
    Comm ** Comm::array;

    Comm::~Comm() {
        if(!rank)
            delete[] array;
    }

    void Comm::send(int to, int size, double * buff, int tag) {
        return;
    }

    void Comm::recv(int from, int size, double * buff, int tag) {
        if(tag>=0) {
            double * src = array[from]->block->edges[tag];
            memcpy(buff, src, size*sizeof(double));
        }

    }

    void Comm::exchange() {
        int sizes[3] = {block->nz*block->ny, block->nz*block->nx, block->nx*block->ny};
        for (int i = 0; i <= 2; ++i) {
            int size = sizes[i];
            // our edge[i]
            double * buff = block->edges[i];

            // edge of our neigbour
            // our left is right for neigbour -> thus 5-i
            double * src = array[neigh[i]]->block->edges[5-i];


            double * temp = new double[size];
            memcpy(temp, buff, size*sizeof(double));
            memcpy(buff, src, size*sizeof(double));
            memcpy(src, temp, size*sizeof(double));
            delete[] temp;
        }
     
    }
#endif
