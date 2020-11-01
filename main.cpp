#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "analytic.h"

#define L_X 1.0
#define L_Y 1.0
#define L_Z 1.0
#define T 0.5

#define N 128
#define K 20

#define PERIOD_X false
#define PERIOD_Y false
#define PERIOD_Z true

const double H_X = L_X / N;
const double H_Y = L_Y / N;
const double H_Z = L_Z / N;

const double TAU = T / K;

const double C_X = TAU / H_X;
const double C_Y = TAU / H_Y;
const double C_Z = TAU / H_Z;

//assumes that layer is double[N+1][N+1][N+1]
double& get3d(double* layer, int i, int j, int k)
{
    return layer[(N+1)*(N+1)*i + (N+1)*j + k];
}


// Laplas operator approximation 
// INCLUDES TAU MULTIPLICATION!
double delta_h(int i, int j, int k, double* curr)
{
    double d_x, d_y, d_z;
    if(PERIOD_X && (i==0 || i==N))
        d_x = get3d(curr, 1, j, k) - 2*get3d(curr, 0, j, k) + get3d(curr, N -1, j, k);
    else
        d_x = get3d(curr, i-1, j, k) - 2*get3d(curr, i, j, k) + get3d(curr, i+1, j, k);

    if(PERIOD_Y && (j==0 || j==N))
        d_y = get3d(curr, i, 1, k) - 2*get3d(curr, i, 0, k) + get3d(curr, i, N-1, k);
    else
        d_y = get3d(curr, i, j-1, k) - 2*get3d(curr, i, j, k) + get3d(curr, i, j+1, k);

    if(PERIOD_Z && (k==0 || k==N))
        d_z = get3d(curr, i, j, 1) - 2*get3d(curr, i, j, 0) + get3d(curr, i, j, N-1);
    else
        d_z = get3d(curr, i, j, k-1) - 2*get3d(curr, i, j, k) + get3d(curr, i, j, k+1);

    return d_x*C_X*C_X + d_y*C_Y*C_Y*0 + d_z*C_Z*C_Z*0;
}

double approx_first(int i, int j, int k, double* zeroth)
{

    if(!PERIOD_X && (i==0 || i==N))
        return 0;

    if(!PERIOD_Y && (j==0 || j==N))
        return 0;

    if(!PERIOD_Z && (k==0 || k==N))
        return 0;


    return get3d(zeroth, i, j, k) + delta_h(i,j,k,zeroth)/2;
}

double approx_next(int i, int j, int k, double* prev, double* curr)
{
    if(!PERIOD_X && (i==0 || i==N))
        return 0;

    if(!PERIOD_Y && (j==0 || j==N))
        return 0;

    if(!PERIOD_Z && (k==0 || k==N))
        return 0;

    if(i ==  N / 2 && j == N / 2 && k == N / 4)
        printf("%9.6f  %9.6f ", get3d(curr, i, j, k) , delta_h(i,j,k,curr));

    return delta_h(i,j,k,curr) + 2*get3d(curr, i, j, k) - get3d(prev, i, j, k);
}

void init_zeroth(double* zeroth)
{
    #pragma omp parallel for
    for(int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j)
            for (int k = 0; k <= N; ++k)
                get3d(zeroth, i, j, k) = phi(L_X, L_Y, L_Z, H_X*i, H_Y*j, H_Z*k);
}

void init_first(double* zeroth, double* first)
{
    #pragma omp parallel for
    for(int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j)
            for (int k = 0; k <= N; ++k)
                // get3d(first, i, j, k) = u_analytical(L_X, L_Y, L_Z, H_X*i, H_Y*j, H_Z*k, TAU);
                get3d(first, i, j, k) = approx_first(i, j, k, zeroth);
}

void calc_next(double* prev, double* curr, double* next)
{
    #pragma omp parallel for
    for(int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j)
            for (int k = 0; k <= N; ++k)
                get3d(next, i, j, k) = approx_next(i, j, k, prev, curr);   
}


double get_error(double* layer, int t)
{
    double max_err=0, temp;
    int mi=0,mj=0,mk=0;
    for(int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j)
            for (int k = 0; k <= N; ++k)
            {
                temp = abs(u_analytical(L_X,L_Y,L_Z, H_X*i, H_Y*j, H_Z*k, t*TAU) - get3d(layer, i, j, k));
                if(temp > max_err)
                {
                    max_err = temp;
                    mi=i;
                    mj=j;
                    mk=k;
                }
            }

    printf("%d %d %d | ", mi, mj, mk);
    return max_err;
}

void print_layer(double* layer, int t)
{
    for(int i = 0; i <= N; ++i)
    {
        for (int j = 0; j <= N; ++j)
        {
            for (int k = 0; k <= N; ++k)
                printf("%7.3f", get3d(layer, i, j, k));
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

void loop()
{
    double * prev = new double[(N+1)*(N+1)*(N+1)];
    double * curr = new double[(N+1)*(N+1)*(N+1)];
    double * next = new double[(N+1)*(N+1)*(N+1)];

    init_zeroth(prev);
    init_first(prev, curr);

    printf("N=%d K=%d H=%.6f TAU=%.6f\n", N, K, H_X, TAU);
    printf("C_X=%.6f C_Y=%.6f C_Z=%.6f\n", C_X, C_Y, C_Z);
    printf("Error #%3d: %8.4f\n", 0, get_error(prev,0));
    printf("Error #%3d: %8.4f\n", 1, get_error(curr,1));


    for(int t=2; t<K; t++)
    {
        printf("%9.6f ", u_analytical(L_X,L_Y,L_Z, H_X*(N/2), H_Y*(N/2), H_Z*(N/4), t*TAU));
        calc_next(prev, curr, next);
        printf("Error #%3d: %8.4f\n", t, get_error(curr,t));

        double * temp = prev;
        prev = curr;
        curr = next;
        next = temp;
    }

    calc_next(prev, curr, next);
    // print_layer(next, K);
    printf("Error #%3d: %8.4f\n", K, get_error(curr,K));

    delete[] prev;
    delete[] curr;
    delete[] next;
}

int main(int argc, char const *argv[])
{
    double start = omp_get_wtime();
    loop();
    printf("%f s elapsed\n", (omp_get_wtime() - start));

    return 0;
}