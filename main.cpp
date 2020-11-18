#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "problem.h"

#define L_X 1.0
#define L_Y 1.0
#define L_Z 1.0
#define T 0.2

#define N 32
#define K 32

const double H_X = L_X / N;
const double H_Y = L_Y / N;
const double H_Z = L_Z / N;

const double TAU = T / K;

const double C_X = TAU / H_X;
const double C_Y = TAU / H_Y;
const double C_Z = TAU / H_Z;

#define BX 1
#define BY 1
#define BZ 1

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
    return (i % n) + (n * (i < 0));
}

struct Block {
    int t = 0;

    double * prev;
    double * curr;
    double * next;

    double * edge_xp;
    double * edge_xm;
    double * edge_yp;
    double * edge_ym;
    double * edge_zp;
    double * edge_zm;

    int sx, ex, nx;
    int sy, ey, ny;
    int sz, ez, nz;

    int px, mx;
    int py, my;
    int pz, mz;

    Block(int rank) {
        int k = rank % BZ;
        int j = rank / BZ % BY;
        int i = rank / BZ / BY;

        sx = (int) round(1.0*N/BX*i);
        ex = (int) round(1.0*N/BX*(i+1));
        nx = ex - sx + 1;
        
        sy = (int) round(1.0*N/BY*j);
        ey = (int) round(1.0*N/BY*(j+1));
        ny = ey - sy + 1;
        
        sz = (int) round(1.0*N/BZ*k);
        ez = (int) round(1.0*N/BZ*(k+1));
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

        edge_xm = new double[ny*nz];
        edge_xp = new double[ny*nz];

        edge_ym = new double[nx*nz];
        edge_yp = new double[nx*nz];

        edge_zm = new double[ny*nx];
        edge_zp = new double[ny*nx];
        printf("sx %d, ex %d\n", sx, ex);
        printf("nx %d, ny %d nz %d\n", nx, ny, nz);
    }

    ~Block() {
        delete[] prev;
        delete[] curr;
        delete[] next;

        delete[] edge_xm;
        delete[] edge_xp;
        delete[] edge_ym;
        delete[] edge_yp;
        delete[] edge_zm;
        delete[] edge_zp;
    }

    void swap() {
        double * temp = prev;
        prev = curr;
        curr = next;
        next = temp;
    }

    double& get(double * layer, int i, int j, int k) {
        if(i==sx-1)
            return edge_xm[nz*(j - sy) + (k - sz)];

        if(i==ex+1)
            return edge_xp[nz*(j - sy) + (k - sz)];

        if(j==sy-1)
            return edge_ym[nz*(i - sx) + (k - sz)];

        if(j==ey+1)
            return edge_yp[nz*(i - sx) + (k - sz)];

        if(k==sz-1)
            return edge_zm[ny*(i - sx) + (j - sy)];

        if(k==ez+1)
            return edge_zp[ny*(i - sx) + (j - sy)];


        return layer[nz*ny*(i - sx) + nz*(j - sy) + (k - sz)];
    }

    void copyAxes(int x, int y, int z, double * from, double * to) {
        int c = 0;
        for(int i = (x<0 ? sx : x); i <= (x<0 ? ex : x); i++)
            for(int j = (y<0 ? sy : y); j <= (y<0 ? ey : y); j++)
                for(int k = (z<0 ? sz : z); k <= (z<0 ? ez : z); k++)                
                    to[c++] = get(from, i,j,k);
    }

    void exchange() {
        copyAxes(sx+1, -1, -1, curr, edge_xp);
        copyAxes(ex-1, -1, -1, curr, edge_xm);
        copyAxes(-1, sy+1, -1, curr, edge_yp);
        copyAxes(-1, ey-1, -1, curr, edge_ym);
        copyAxes(-1, -1, sz+1, curr, edge_zp);
        copyAxes(-1, -1, ez-1, curr, edge_zm);
    }

    void init0() {
        #pragma omp parallel for
        for(int i = sx; i <= ex; i++)
            for(int j = sy; j <= ey;  j++)
                for(int k = sz; k <= ez; k++)
                    get(next, i, j, k) = phi(L_X, L_Y, L_Z, H_X*i, H_Y*j, H_Z*k);

    }

    double delta(int i, int j, int k, double* curr) {
        double d_x, d_y, d_z;

        d_x = (get(curr, i+1, j, k) - get(curr, i, j, k)) * C_X + (get(curr, i-1, j, k) - get(curr, i, j, k)) * C_X;
        d_y = (get(curr, i, j-1, k) - get(curr, i, j, k)) * C_Y + (get(curr, i, j+1, k) - get(curr, i, j, k)) * C_Y;
        d_z = (get(curr, i, j, k-1) - get(curr, i, j, k)) * C_Z + (get(curr, i, j, k+1) - get(curr, i, j, k)) * C_Z;

        return d_x*C_X + d_y*C_Y + d_z*C_Z;
    }

    void init1() {
        swap();
        exchange();
        #pragma omp parallel for
        for(int i = sx; i <= ex; i++)
            for(int j = sy; j <= ey; j++)
                for(int k = sz; k <= ez; k++) {
                    if(onConstEdge(i,j,k))
                        get(next, i, j, k) = 0;
                    else
                        get(next, i, j, k) = get(curr, i, j, k) + delta(i,j,k,curr)/2;
                }
        t++;        
    }

    void calcNext() {
        swap();
        exchange();
        #pragma omp parallel for
        for(int i = sx; i <= ex; i++)
            for(int j = sy; j <= ey; j++)
                for(int k = sz; k <= ez; k++) {
                    if(onConstEdge(i,j,k))
                        get(next, i, j, k) = 0;
                    else
                        get(next, i, j, k) = get(curr, i, j, k) + (delta(i,j,k,curr) + get(curr, i, j, k) - get(prev, i, j, k));

                }
        t++;
    }

    double get_error()
    {
        double max_err=0, temp;

        for(int i = 0; i <= N; ++i)
            for (int j = 0; j <= N; ++j)
                for (int k = 0; k <= N; ++k)
                {
                    temp = std::abs(u_analytical(L_X,L_Y,L_Z, H_X*i, H_Y*j, H_Z*k, t*TAU) - get(next, i, j, k));
                    if(temp > max_err)
                        max_err = temp;
                }

        return max_err;
    }

    void print_layer()
    {
        for(int i = 0; i <= N; ++i)
        {
            for (int j = 0; j <= N; ++j)
            {
                for (int k = 0; k <= N; ++k)
                    printf("%7.3f", get(next, i, j, k));
                printf("\n");
                printf("\033[1;30m");
                for (int k = 0; k <= N; ++k)
                    printf("%7.3f", u_analytical(L_X,L_Y,L_Z, H_X*i, H_Y*j, H_Z*k, t*TAU));
                printf(" ***\n");
                printf("\033[0m");
            }
            printf("\n\n");
        }
    }
};

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
