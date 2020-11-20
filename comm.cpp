#include "comm.hpp"
#include "block.hpp"
#include <string.h>

// define static member
Comm ** Comm::array;

Comm::Comm(int rank, int size, Block * block) {
    this->rank = rank;
    this->size = size;
    this->block = block;

    if(!rank)
        array = new Comm*[size];

    array[rank] = this;
}

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

void Comm::swap(int with, int size, double * buff, int tag) {
   if(tag >= 0 && tag <= 2) {
        double * src = array[with]->block->edges[tag];
        double * temp = new double[size];
        memcpy(temp, buff, size*sizeof(double));
        memcpy(buff, src, size*sizeof(double));
        memcpy(src, temp, size*sizeof(double));
        delete[] temp;
    }
 
}