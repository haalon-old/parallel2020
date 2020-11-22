#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "block.hpp"
#include "comm.hpp"
#include "problem.h"

char onConstEdge(int i, int j, int k) {
    if(!PERIOD_X && (i==0 || i==N))
        return 1;

    if(!PERIOD_Y && (j==0 || j==N))
        return 1;

    if(!PERIOD_Z && (k==0 || k==N))
        return 1;

    return 0;
}

int mod(int i, int n) {
    return ((i % n) + n) % n;
}

Block::Block() {}

Block::Block(int rank) {
    int k = rank % BZ;
    int j = rank / BZ % BY;
    int i = rank / BZ / BY;

    // if axis is periodic  
    // the first and the last node on that axis
    // is equivalent

    sx = (int) round(1.0*N/BX*i);
    ex = (int) round(1.0*N/BX*(i+1))-1;
    nx = ex - sx + 1;
    
    sy = (int) round(1.0*N/BY*j);
    ey = (int) round(1.0*N/BY*(j+1))-1;
    ny = ey - sy + 1;
    
    sz = (int) round(1.0*N/BZ*k);
    ez = (int) round(1.0*N/BZ*(k+1))-1;
    nz = ez - sz + 1;

    // ranks of adjacent blocks
    pz = mod(k+1, BZ) + BZ*j + BZ*BY*i;
    mz = mod(k-1, BZ) + BZ*j + BZ*BY*i;

    py = k + BZ*mod(j+1, BY) + BZ*BY*i;
    my = k + BZ*mod(j-1, BY) + BZ*BY*i;

    px = k + BZ*j + BZ*BY*mod(i+1, BX);
    mx = k + BZ*j + BZ*BY*mod(i-1, BX);

    prev = new double[nx*ny*nz];
    curr = new double[nx*ny*nz];
    next = new double[nx*ny*nz];

    // sum of indexes of opposite sides is 5
    edges[0] = new double[ny*nz]; // x-
    edges[5] = new double[ny*nz]; // x+

    edges[1] = new double[nx*nz]; // y-
    edges[4] = new double[nx*nz]; // y+

    edges[2] = new double[ny*nx]; // z-
    edges[3] = new double[ny*nx]; // z+

    printf("\tsx %d, ex %d\n", sx, ex);
    printf("\tsy %d, ey %d\n", sy, ey);
    printf("\tsz %d, ez %d\n", sz, ez);
    printf("\tnx %d, ny %d nz %d\n", nx, ny, nz);

    printf("\tx %d %d, y %d %d, z %d %d\n", mx, px, my, py, mz, pz);
}

Block::~Block() {
    delete[] prev;
    delete[] curr;
    delete[] next;

    for(int i=0; i<6; i++)
        delete[] edges[i];
}

void Block::swap() {
    double * temp = prev;
    prev = curr;
    curr = next;
    next = temp;
}

double& Block::get(double * layer, int i, int j, int k) {
    if(i==sx-1)
        return edges[0][nz*(j - sy) + (k - sz)];

    if(i==ex+1)
        return edges[5][nz*(j - sy) + (k - sz)];

    if(j==sy-1)
        return edges[1][nz*(i - sx) + (k - sz)];

    if(j==ey+1)
        return edges[4][nz*(i - sx) + (k - sz)];

    if(k==sz-1)
        return edges[2][ny*(i - sx) + (j - sy)];

    if(k==ez+1)
        return edges[3][ny*(i - sx) + (j - sy)];


    return layer[nz*ny*(i - sx) + nz*(j - sy) + (k - sz)];
}

void Block::copyAxes(int x, int y, int z, double * from, double * to) {
    int c = 0;
    for(int i = (x<0 ? sx : x); i <= (x<0 ? ex : x); i++)
        for(int j = (y<0 ? sy : y); j <= (y<0 ? ey : y); j++)
            for(int k = (z<0 ? sz : z); k <= (z<0 ? ez : z); k++)                
                to[c++] = get(from, i,j,k);
}

void Block::prepare() {
    // 
    copyAxes(ex, -1, -1, next, edges[5]);
    copyAxes(sx, -1, -1, next, edges[0]);
    
    copyAxes(-1, ey, -1, next, edges[4]);
    copyAxes(-1, sy, -1, next, edges[1]);
    
    copyAxes(-1, -1, ez, next, edges[3]);
    copyAxes(-1, -1, sz, next, edges[2]);
}

void Block::exchange(Comm * comm) {

    // we send edge[5], we expect edge[0] from the x+ neighbour
    comm->swap(px, ny*nz, edges[5], 0);
    comm->swap(mx, ny*nz, edges[0], 5);

    comm->swap(py, nx*nz, edges[4], 1);
    comm->swap(my, nx*nz, edges[1], 4);

    comm->swap(pz, nx*ny, edges[3], 2);
    comm->swap(pz, nx*ny, edges[2], 3);
}

void Block::init0() {
    // #pragma omp parallel for
    for(int i = sx; i <= ex; i++)
        for(int j = sy; j <= ey;  j++)
            for(int k = sz; k <= ez; k++)
                get(next, i, j, k) = phi(L_X, L_Y, L_Z, H_X*i, H_Y*j, H_Z*k);
    prepare();
}

double Block::delta(int i, int j, int k, double* curr) {
    double d_x, d_y, d_z;

    d_x = (get(curr, i+1, j, k) - get(curr, i, j, k)) * C_X + (get(curr, i-1, j, k) - get(curr, i, j, k)) * C_X;
    d_y = (get(curr, i, j-1, k) - get(curr, i, j, k)) * C_Y + (get(curr, i, j+1, k) - get(curr, i, j, k)) * C_Y;
    d_z = (get(curr, i, j, k-1) - get(curr, i, j, k)) * C_Z + (get(curr, i, j, k+1) - get(curr, i, j, k)) * C_Z;

    return d_x*C_X + d_y*C_Y + d_z*C_Z;
}

void Block::init1() {
    swap();
    // exchange(comm);
    // #pragma omp parallel for
    for(int i = sx; i <= ex; i++)
        for(int j = sy; j <= ey; j++)
            for(int k = sz; k <= ez; k++) {
                if(onConstEdge(i,j,k))
                    get(next, i, j, k) = 0;
                else
                    get(next, i, j, k) = get(curr, i, j, k) + delta(i,j,k,curr)/2;
            }
    t++;
    prepare();
}

void Block::calcNext() {
    swap();
    // exchange(comm);
    // #pragma omp parallel for
    for(int i = sx; i <= ex; i++)
        for(int j = sy; j <= ey; j++)
            for(int k = sz; k <= ez; k++) {
                if(onConstEdge(i,j,k))
                    get(next, i, j, k) = 0;
                else
                    get(next, i, j, k) = get(curr, i, j, k) + (delta(i,j,k,curr) + get(curr, i, j, k) - get(prev, i, j, k));

            }
    t++;
    prepare();
}

double Block::get_error()
{
    double max_err=0, temp;

    for(int i = sx; i <= ex; ++i)
        for (int j = sy; j <= ey; ++j)
            for (int k = sz; k <= ez; ++k)
            {
                temp = std::abs(u_analytical(L_X,L_Y,L_Z, H_X*i, H_Y*j, H_Z*k, t*TAU) - get(next, i, j, k));
                if(temp > max_err)
                    max_err = temp;
            }

    return max_err;
}

void Block::print_layer()
{
    for(int i = sx; i <= ex; ++i)
    {
        for (int j = sy; j <= ey; ++j)
        {
            for (int k = sz; k <= ez; ++k)
                printf("%7.3f", get(next, i, j, k));
            printf("\n");
            printf("\033[1;30m");
            for (int k = sz; k <= ez; ++k)
                printf("%7.3f", u_analytical(L_X,L_Y,L_Z, H_X*i, H_Y*j, H_Z*k, t*TAU));
            printf(" ***\n");
            printf("\033[0m");
        }
        printf("\n\n");
    }
}