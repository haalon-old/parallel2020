#pragma once

struct Comm;

struct Block {
    int t = 0;

    double * prev;
    double * curr;
    double * next;

    double * edges[6];

    int sx, ex, nx;
    int sy, ey, ny;
    int sz, ez, nz;

    int px, mx;
    int py, my;
    int pz, mz;

    Block();
    Block(int rank);
    ~Block();

    double& get(double * layer, int i, int j, int k);

    void swap();

    void copyAxes(int x, int y, int z, double * from, double * to);
    void prepare();
    void exchange(Comm * comm);
    
    double delta(int i, int j, int k, double* curr);

    void init0();
    void init1();
    void calcNext();

    double get_error();
    void print_layer();
};