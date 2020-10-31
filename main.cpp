#include <stdio.h>
#include <math.h>
#include "analytic.h"

#define L_X 10.0
#define L_Y 10.0
#define L_Z 10.0
#define T 1.0

#define N 128
#define K 128

#define H_X L_X / N
#define H_Y L_Y / N
#define H_Z L_Z / N

#define TAU T / K

#define PERIOD_X false
#define PERIOD_Y false
#define PERIOD_Z true


double& get3d(double* layer, int i, int j, int k)
{
    return layer[(N+1)*(N+1)*i + (N+1)*j + k];
}


// Laplas operator approximation
double delta_h(int i, int j, int k, double* curr)
{
    double d_x, d_y, d_z;
    if(PERIOD_X && (i==0 || i==N))
        d_x = (get3d(curr, 1, j, k) - 2*get3d(curr, 0, j, k) + get3d(curr, N -1, j, k)) / (H_X*H_X);
    else
        d_x = (get3d(curr, i-1, j, k) - 2*get3d(curr, i, j, k) + get3d(curr, i+1, j, k)) / (H_X*H_X);

    if(PERIOD_Y && (j==0 || j==N))
        d_y = (get3d(curr, i, 1, k) - 2*get3d(curr, i, 0, k) + get3d(curr, i, N-1, k)) / (H_Y*H_Y);
    else
        d_y = (get3d(curr, i, j-1, k) - 2*get3d(curr, i, j, k) + get3d(curr, i, j+1, k)) / (H_Y*H_Y);

    if(PERIOD_Z && (k==0 || k==N))
        d_z = (get3d(curr, i, j, 1) - 2*get3d(curr, i, j, 0) + get3d(curr, i, j, N-1)) / (H_Z*H_Z);
    else
        d_z = (get3d(curr, i, j, k-1) - 2*get3d(curr, i, j, k) + get3d(curr, i, j, k+1)) / (H_Z*H_Z);

    return d_x + d_y + d_z;
}

double approx_first(int i, int j, int k, double* zeroth)
{
#if !PERIOD_X
    if(i==0 || i==N)
        return 0;
#endif

#if !PERIOD_Y
    if(j==0 || j==N)
        return 0;
#endif

#if !PERIOD_Z
    if(k==0 || k==N)
        return 0;
#endif

    return get3d(zeroth, i, j, k) + TAU*TAU/2*delta_h(i,j,k,zeroth);
}

double approx_next(int i, int j, int k, double* prev, double* curr)
{
#if !PERIOD_X
    if(i==0 || i==N)
        return 0;
#endif

#if !PERIOD_Y
    if(j==0 || j==N)
        return 0;
#endif

#if !PERIOD_Z
    if(k==0 || k==N)
        return 0;
#endif

    return TAU*TAU*delta_h(i,j,k,curr) + 2*get3d(curr, i, j, k) - get3d(prev, i, j, k);
}

void init_zeroth(double* zeroth)
{
    for(int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j)
            for (int k = 0; k <= N; ++k)
                get3d(zeroth, i, j, k) = phi(L_X, L_Y, L_Z, H_X*i, H_Y*j, H_Z*k);
}

void init_first(double* zeroth, double* first)
{
    for(int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j)
            for (int k = 0; k <= N; ++k)
                get3d(first, i, j, k) = approx_first(i, j, k, zeroth);
}

void calc_next(double* prev, double* curr, double* next)
{
    for(int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j)
            for (int k = 0; k <= N; ++k)
                get3d(next, i, j, k) = approx_next(i, j, k, prev, curr);   
}


double get_error(double* layer, int t)
{
    double max_err=0, temp;
    for(int i = 0; i <= N; ++i)
        for (int j = 0; j <= N; ++j)
            for (int k = 0; k <= N; ++k)
            {
                temp = abs(u_analytical(L_X,L_Y,L_Z, H_X*i, H_Y*j, H_Z*k, t*TAU) - get3d(layer, i, j, k));
                max_err = max_err > temp ? max_err : temp;
            }

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

            for (int k = 0; k <= N; ++k)
                printf("%7.3f", u_analytical(L_X,L_Y,L_Z, H_X*i, H_Y*j, H_Z*k, t*TAU));
            printf(" ***\n");
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

    printf("N=%d K=%d\n", N, K);
    printf("Error #0: %8.4f\n", get_error(prev,0));
    printf("Error #1: %8.4f\n", get_error(curr,1));


    for(int t=2; t<K; t++)
    {
        calc_next(prev, curr, next);
        printf("Error #%d: %8.4f\n", t, get_error(curr,t));

        double * temp = prev;
        prev = curr;
        curr = next;
        next = temp;
    }

    calc_next(prev, curr, next);
    // print_layer(next, K);
    printf("Error #%d: %8.4f\n", K, get_error(curr,K));

    delete[] prev;
    delete[] curr;
    delete[] next;
}

int main(int argc, char const *argv[])
{
    loop();
    return 0;
}