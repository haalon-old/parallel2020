#pragma once

struct Block;

struct Comm
{
    int rank, size;
    static Comm ** array;
    Block * block;

    Comm();
    Comm(int rank, int size, Block * block);
    ~Comm();

    void send(int to, int size, double * buff, int tag);    
    void recv(int from, int size, double * buff, int tag);

    void swap(int with, int size, double * buff, int tag);
};