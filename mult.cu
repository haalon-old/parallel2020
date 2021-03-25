#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include <thrust/extrema.h>
#include <thrust/execution_policy.h>

#include "problem.h"
#include "block.hpp"

#define REDUCTION_BLOCK_SIZE 1024
#define THREADSINBLOCK 512
#define THREAD_AXIS 8

__constant__ int N_D;

__constant__ double C_X_D;
__constant__ double C_Y_D;
__constant__ double C_Z_D;

__constant__ int sx_d, ex_d, nx_d;
__constant__ int sy_d, ey_d, ny_d;
__constant__ int sz_d, ez_d, nz_d;

__constant__ double H_X_D, H_Y_D, H_Z_D, TAU_D;
__constant__ double L_X_D, L_Y_D, L_Z_D;

__constant__ double * prev_d;
__constant__ double * curr_d;
__constant__ double * next_d;

__constant__ double * edges_d[6];
__constant__ double * new_edges_d[6];

size_t edge_sizes[6];

__device__ double u_analytical_d(double l_x, double l_y, double l_z, double x, double y, double z, double t)
{
    double a_t = PI * sqrt(1/(l_x*l_x) + 4/(l_y*l_y) + 9/(l_z*l_z));
    return sin(PI*x/l_x) * sin(2*PI*y/l_y) * sin(3*PI*z/l_z) * cos(a_t*t);
}

#define SAFE_CALL( CallInstruction ) { \
    cudaError_t cuerr = CallInstruction; \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error: %s at call \"" #CallInstruction "\"\n", cudaGetErrorString(cuerr)); \
        throw "error in CUDA API function, aborting..."; \
    } \
}


#define SAFE_KERNEL_CALL( KernelCallInstruction ){ \
    KernelCallInstruction; \
    cudaError_t cuerr = cudaGetLastError(); \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error in kernel launch: %s at kernel \"" #KernelCallInstruction "\"\n", cudaGetErrorString(cuerr)); \
        throw "error in CUDA kernel launch, aborting..."; \
    } \
    cuerr = cudaDeviceSynchronize(); \
    if(cuerr != cudaSuccess) { \
        printf("CUDA error in kernel execution: %s at kernel \"" #KernelCallInstruction "\"\n", cudaGetErrorString(cuerr)); \
        throw "error in CUDA kernel execution, aborting..."; \
    } \
}

extern void initDevice(Block * b){
    // set extern consts
    SAFE_CALL(cudaMemcpyToSymbol(N_D, &N, sizeof(int)));

    SAFE_CALL(cudaMemcpyToSymbol(C_X_D, &C_X, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(C_Y_D, &C_Y, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(C_Z_D, &C_Z, sizeof(double)));

    SAFE_CALL(cudaMemcpyToSymbol(H_X_D, &H_X, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(H_Y_D, &H_Y, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(H_Z_D, &H_Z, sizeof(double)));

    SAFE_CALL(cudaMemcpyToSymbol(L_X_D, &L_X, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(L_Y_D, &L_Y, sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(L_Z_D, &L_Z, sizeof(double)));

    SAFE_CALL(cudaMemcpyToSymbol(TAU_D, &TAU, sizeof(double)));

    // set block parameters
    SAFE_CALL(cudaMemcpyToSymbol(sx_d, &b->sx, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(ex_d, &b->ex, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(nx_d, &b->nx, sizeof(int)));

    SAFE_CALL(cudaMemcpyToSymbol(sy_d, &b->sy, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(ey_d, &b->ey, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(ny_d, &b->ny, sizeof(int)));

    SAFE_CALL(cudaMemcpyToSymbol(sz_d, &b->sz, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(ez_d, &b->ez, sizeof(int)));
    SAFE_CALL(cudaMemcpyToSymbol(nz_d, &b->nz, sizeof(int)));


    // set arrays
    int size = b->nx * b->ny * b->nz;
    int sizeX = b->ny * b->nz, sizeY = b->nx * b->nz, sizeZ = b->nx * b->ny;
    edge_sizes[0] = sizeX; edge_sizes[1] = sizeY; edge_sizes[2] = sizeZ;
    edge_sizes[3] = sizeZ; edge_sizes[4] = sizeY; edge_sizes[5] = sizeX;
    double * temp_host;

    // 1) allocate memory, temp_host will have pointer to dev memory
    // (not needed anymore!) 2) copy host data to device data 
    // 3) copy the pointer itself to the device
    SAFE_CALL(cudaMalloc((void**)&temp_host, size * sizeof(double)));       
    SAFE_CALL(cudaMemcpyToSymbol(prev_d, &temp_host, sizeof(double *)));

    SAFE_CALL(cudaMalloc((void**)&temp_host, size * sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(curr_d, &temp_host, sizeof(double *)));

    SAFE_CALL(cudaMalloc((void**)&temp_host, size * sizeof(double)));
    SAFE_CALL(cudaMemcpyToSymbol(next_d, &temp_host, sizeof(double *)));

    double * host_edges_temp[6];

    for (int i = 0; i < 6; ++i) 
        SAFE_CALL(cudaMalloc((void**)&host_edges_temp[i], edge_sizes[i] * sizeof(double)));
    
    SAFE_CALL(cudaMemcpyToSymbol(edges_d, &host_edges_temp, 6 * sizeof(double *)));

    for (int i = 0; i < 6; ++i) 
        SAFE_CALL(cudaMalloc((void**)&host_edges_temp[i], edge_sizes[i] * sizeof(double)));
    
    SAFE_CALL(cudaMemcpyToSymbol(new_edges_d, &host_edges_temp, 6 * sizeof(double *)));
}

void swap() {
    double * next, curr, prev;
    SAFE_CALL(cudaMemcpyFromSymbol(&next, next_d, sizeof(double *)));
    SAFE_CALL(cudaMemcpyFromSymbol(&curr, curr_d, sizeof(double *)));
    SAFE_CALL(cudaMemcpyFromSymbol(&prev, prev_d, sizeof(double *)));

    SAFE_CALL(cudaMemcpyToSymbol(curr_d, &next, sizeof(double *)));
    SAFE_CALL(cudaMemcpyToSymbol(prev_d, &curr, sizeof(double *)));
    SAFE_CALL(cudaMemcpyToSymbol(next_d, &prev, sizeof(double *)));
}

__device__ char onConstEdge_D(int i, int j, int k) {
    if(!PERIOD_X && (i==0 || i==N_D))
        return 1;

    if(!PERIOD_Y && (j==0 || j==N_D))
        return 1;

    if(!PERIOD_Z && (k==0 || k==N_D))
        return 1;

    return 0;
}

__device__ double& get(double * layer,  int i, int j, int k)
{
    if(i==sx_d-1)
        return edges_d[0][nz_d*(j - sy_d) + (k - sz_d)];

    if(i==ex_d+1)
        return edges_d[5][nz_d*(j - sy_d) + (k - sz_d)];

    if(j==sy_d-1)
        return edges_d[1][nz_d*(i - sx_d) + (k - sz_d)];

    if(j==ey_d+1)
        return edges_d[4][nz_d*(i - sx_d) + (k - sz_d)];

    if(k==sz_d-1)
        return edges_d[2][ny_d*(i - sx_d) + (j - sy_d)];

    if(k==ez_d+1)
        return edges_d[3][ny_d*(i - sx_d) + (j - sy_d)];


    return layer[nz_d*ny_d*(i - sx_d) + nz_d*(j - sy_d) + (k - sz_d)];
}

__device__ void setNewEdges(double val, int i, int j, int k) {
    if(i==sx_d)
        new_edges_d[0][nz_d*(j - sy_d) + (k - sz_d)] = val;

    if(i==ex_d)
        new_edges_d[5][nz_d*(j - sy_d) + (k - sz_d)] = val;

    if(j==sy_d)
        new_edges_d[1][nz_d*(i - sx_d) + (k - sz_d)] = val;

    if(j==ey_d)
        new_edges_d[4][nz_d*(i - sx_d) + (k - sz_d)] = val;

    if(k==sz_d)
        new_edges_d[2][ny_d*(i - sx_d) + (j - sy_d)] = val;

    if(k==ez_d)
        new_edges_d[3][ny_d*(i - sx_d) + (j - sy_d)] = val;
}

__device__ double delta(double * curr, int i, int j, int k) {
    double d_x, d_y, d_z;

    d_x = (get(curr, i+1, j, k) - get(curr, i, j, k)) * C_X_D + (get(curr, i-1, j, k) - get(curr, i, j, k)) * C_X_D;
    d_y = (get(curr, i, j-1, k) - get(curr, i, j, k)) * C_Y_D + (get(curr, i, j+1, k) - get(curr, i, j, k)) * C_Y_D;
    d_z = (get(curr, i, j, k-1) - get(curr, i, j, k)) * C_Z_D + (get(curr, i, j, k+1) - get(curr, i, j, k)) * C_Z_D;

    return d_x*C_X_D + d_y*C_Y_D + d_z*C_Z_D;
}



__global__ void __calc_n__(int num) {
    int gridOffsetX = blockDim.x * blockIdx.x;
    int gridOffsetY = blockDim.y * blockIdx.y;
    int gridOffsetZ = blockDim.z * blockIdx.z;

    // shift by s(.)_d, so i j k are indexing from the start of the block
    int i = threadIdx.x + gridOffsetX + sx_d;
    int j = threadIdx.y + gridOffsetY + sy_d;
    int k = threadIdx.z + gridOffsetZ + sz_d;

    // it shouldnt really happen, unless we use weird grid sizes
    if(i > ex_d || j > ey_d || k > ez_d)
        return;

    double val;

    switch(num) {
        case  0: val = u_analytical_d(L_X_D, L_Y_D, L_Z_D, H_X_D*i, H_Y_D*j, H_Z_D*k, 0); break;
        case  1: val = onConstEdge_D(i,j,k) ? 0 : get(curr_d, i, j, k) + delta(curr_d, i,j,k)/2.0; break;
        default: val = onConstEdge_D(i,j,k) ? 0 : 2*get(curr_d, i, j, k) + delta(curr_d, i, j, k) - get(prev_d, i, j, k); break;
    }

    get(next_d, i,j,k) = val;
    setNewEdges(val, i,j,k);

}

extern void launch_calc(Block * b, int num)
{
    double * temp_edges[6];

    // don't needed on the zeroth step
    if(num) {
        // same place as in block
        swap();

        // copy edges we recieved to the gpu
        SAFE_CALL(cudaMemcpyFromSymbol(&temp_edges, edges_d, 6 * sizeof(double *)));
        for (int i = 0; i < 6; ++i) {
            SAFE_CALL(cudaMemcpy((void*)temp_edges[i], (void*)b->edges[i], edge_sizes[i] * sizeof(double), cudaMemcpyHostToDevice));
        }
    }

    dim3 blockDim = dim3(THREAD_AXIS,THREAD_AXIS,THREAD_AXIS);
    int gridDimX = (b->nx - 1)/THREAD_AXIS + 1;
    int gridDimY = (b->ny - 1)/THREAD_AXIS + 1;
    int gridDimZ = (b->nz - 1)/THREAD_AXIS + 1;
    dim3 gridDim = dim3(gridDimX, gridDimY, gridDimZ);

    // get pointer to next array on device, and put it into temp
    SAFE_CALL(cudaMemcpyFromSymbol(&temp_edges, new_edges_d, 6 * sizeof(double *)));

    SAFE_KERNEL_CALL((__calc_n__<<<gridDim, blockDim>>>(num)));    


    // SAFE_CALL(cudaMemcpy((void*)b->next, (void*)temp, size * sizeof(double), cudaMemcpyDeviceToHost));
    for (int i = 0; i < 6; ++i)
        SAFE_CALL(cudaMemcpy((void*)b->edges[i], (void*)temp_edges[i], edge_sizes[i] * sizeof(double), cudaMemcpyDeviceToHost));
}


__device__ double diff(unsigned int indx, int t) {
    unsigned int k = indx % nz_d        + sz_d;
    unsigned int j = indx / nz_d % ny_d + sy_d;
    unsigned int i = indx / nz_d / ny_d + sx_d;
    
    return abs(next_d[indx] - u_analytical_d(L_X_D,L_Y_D,L_Z_D, H_X_D*i, H_Y_D*j, H_Z_D*k, t*TAU_D));
}

__global__ void __err__(int t) {
    int gridOffsetX = blockDim.x * blockIdx.x;
    int gridOffsetY = blockDim.y * blockIdx.y;
    int gridOffsetZ = blockDim.z * blockIdx.z;

    // shift by s(.)_d, so i j k are indexing from the start of the block
    int i = threadIdx.x + gridOffsetX + sx_d;
    int j = threadIdx.y + gridOffsetY + sy_d;
    int k = threadIdx.z + gridOffsetZ + sz_d;

    // it shouldnt really happen, unless we use weird grid sizes
    if(i > ex_d || j > ey_d || k > ez_d)
        return;

    get(prev_d, i,j,k) = abs(get(next_d, i,j,k) - u_analytical_d(L_X_D,L_Y_D,L_Z_D, H_X_D*i, H_Y_D*j, H_Z_D*k, t*TAU_D));

}

extern double launch_err(Block * b) 
{
    // double * dev_out;
    // double out;
    // SAFE_CALL(cudaMalloc((void**)&dev_out, sizeof(double)));

    // SAFE_KERNEL_CALL((__err__<<<1, REDUCTION_BLOCK_SIZE>>>(b->t, dev_out)));

    // SAFE_CALL(cudaMemcpy(&out, dev_out, sizeof(double), cudaMemcpyDeviceToHost));
    // SAFE_CALL(cudaFree(dev_out));

    // return out;

    double * temp;
    double * res;
    double out;
    int size = b->nx * b->ny * b->nz;

    dim3 blockDim = dim3(THREAD_AXIS,THREAD_AXIS,THREAD_AXIS);
    int gridDimX = (b->nx - 1)/THREAD_AXIS + 1;
    int gridDimY = (b->ny - 1)/THREAD_AXIS + 1;
    int gridDimZ = (b->nz - 1)/THREAD_AXIS + 1;
    dim3 gridDim = dim3(gridDimX, gridDimY, gridDimZ);

    SAFE_CALL(cudaMemcpyFromSymbol(&temp, prev_d, sizeof(double *)));
    SAFE_KERNEL_CALL((__err__<<<gridDim, blockDim>>>(b->t)));
    

    res = thrust::max_element(thrust::device, temp, temp + size);
    SAFE_CALL(cudaMemcpy(&out, res, sizeof(double), cudaMemcpyDeviceToHost));

    return out;
}

extern void freeDevice() {
    double * temp;
    double * temp_edges[6];

    SAFE_CALL(cudaMemcpyFromSymbol(&temp, prev_d, sizeof(double *)));
    SAFE_CALL(cudaFree(temp));

    SAFE_CALL(cudaMemcpyFromSymbol(&temp, curr_d, sizeof(double *)));
    SAFE_CALL(cudaFree(temp));

    SAFE_CALL(cudaMemcpyFromSymbol(&temp, next_d, sizeof(double *)));
    SAFE_CALL(cudaFree(temp));

    SAFE_CALL(cudaMemcpyFromSymbol(&temp_edges, edges_d, 6 * sizeof(double *)));
    for(int i=0; i<6; i++)
        SAFE_CALL(cudaFree(temp_edges[i]));
}