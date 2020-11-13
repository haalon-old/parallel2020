#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "analytic.h"

#define L_X 1.0
#define L_Y 1.0
#define L_Z 1.0
#define T 0.2

#define N 128
#define K 128

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

static int t_global=0;

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
        d_x = (get3d(curr, 1, j, k) - get3d(curr, 0, j, k)) * C_X + (get3d(curr, N -1, j, k) - get3d(curr, 0, j, k)) * C_X;
    else
        d_x = (get3d(curr, i+1, j, k) - get3d(curr, i, j, k)) * C_X + (get3d(curr, i-1, j, k) - get3d(curr, i, j, k)) * C_X;

    if(PERIOD_Y && (j==0 || j==N))
        d_y = (get3d(curr, i, 1, k) - get3d(curr, i, 0, k)) * C_Y + (get3d(curr, i, N-1, k) - get3d(curr, i, 0, k)) * C_Y;
    else
        d_y = (get3d(curr, i, j-1, k) - get3d(curr, i, j, k)) * C_Y + (get3d(curr, i, j+1, k) - get3d(curr, i, j, k)) * C_Y;

    if(PERIOD_Z && (k==0 || k==N))
        d_z = (get3d(curr, i, j, 1) - get3d(curr, i, j, 0)) * C_Z + (get3d(curr, i, j, N-1) - get3d(curr, i, j, 0)) * C_Z;
    else
        d_z = (get3d(curr, i, j, k-1) - get3d(curr, i, j, k)) * C_Z + (get3d(curr, i, j, k+1) - get3d(curr, i, j, k)) * C_Z;

    return d_x*C_X + d_y*C_Y + d_z*C_Z;
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

    // if(t_global >= 10 && i==65 && j==11 && k==24)
    // {
    //     double next = delta_h(i,j,k,curr) + 2*get3d(curr, i, j, k) - get3d(prev, i, j, k);
    //     double d_x = (get3d(curr, i+1, j, k) - get3d(curr, i, j, k)) * C_X + (get3d(curr, i-1, j, k) - get3d(curr, i, j, k)) * C_X;
    //     double d_y = (get3d(curr, i, j-1, k) - get3d(curr, i, j, k)) * C_Y + (get3d(curr, i, j+1, k) - get3d(curr, i, j, k)) * C_Y;
    //     double d_z = (get3d(curr, i, j, k-1) - get3d(curr, i, j, k)) * C_Z + (get3d(curr, i, j, k+1) - get3d(curr, i, j, k)) * C_Z;
    //     printf("\n\n%d %d %d (%d)\n",i, j, k, t_global);
    //     printf("prev real (appr) value: %.17f (%.17f)\n", u_analytical(L_X,L_Y,L_Z, H_X*i, H_Y*j, H_Z*k, t_global*TAU - TAU - TAU), get3d(prev, i, j, k));
    //     printf("curr real (appr) value: %.17f (%.17f)\n", u_analytical(L_X,L_Y,L_Z, H_X*i, H_Y*j, H_Z*k, t_global*TAU - TAU), get3d(curr, i, j, k));
    //     printf("next real (appr) value: %.17f (%.17f)\n", u_analytical(L_X,L_Y,L_Z, H_X*i, H_Y*j, H_Z*k, t_global*TAU), next);        
    //     printf("%s: %.17f\n", "xyz", get3d(curr, i, j, k));
    //     printf("%s: %.17f\n", "x+1", get3d(curr, i+1, j, k));
    //     printf("%s: %.17f\n", "x-1", get3d(curr, i-1, j, k));
    //     printf("%s: %.17f\n", "y+1", get3d(curr, i, j+1, k));
    //     printf("%s: %.17f\n", "y-1", get3d(curr, i, j-1, k));
    //     printf("%s: %.17f\n", "z+1", get3d(curr, i, j, k+1));
    //     printf("%s: %.17f\n", "z-1", get3d(curr, i, j, k-1));
    //     printf("%s: %.17f\n", "z-1", get3d(curr, i, j, k-1));
    //     printf("dx %.17f\n", d_x);
    //     printf("dy %.17f\n", d_y);
    //     printf("dz %.17f\n", d_z);
    //     printf("dlt: %.17f\n", delta_h(i,j,k,curr));

    //     // if(delta_h(i,j,k,curr) > 8)
    //     //     abort();
    // }

    return get3d(curr, i, j, k) + (delta_h(i,j,k,curr) + get3d(curr, i, j, k) - get3d(prev, i, j, k));
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
    printf("Error #%3d: %.17f\n", 0, get_error(prev,0));
    printf("Error #%3d: %.17f\n", 1, get_error(curr,1));


    for(int t=2; t<K; t++)
    {
        t_global = t;
        calc_next(prev, curr, next);
        printf("Error #%3d: %.17f\n", t, get_error(next,t));

        double * temp = prev;
        prev = curr;
        curr = next;
        next = temp;
    }

    calc_next(prev, curr, next);
    // print_layer(next, K);
    printf("Error #%3d: %.17f\n", K, get_error(next,K));

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