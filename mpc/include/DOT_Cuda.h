/**
 * @fileoverview Copyright (c) 2019-2022, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#ifdef __MY_CUDA

// Zeta block for coordinates vector
#define BLOCKSIZE 1024

// Taken from:
#define CHECK(call)                                          \
  {                                                          \
    const cudaError_t error = call;                          \
    if (error != cudaSuccess) {                              \
      fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__); \
      fprintf(stderr, "code: %d, reason: %s\n", error,       \
              cudaGetErrorString(error));                    \
      exit(1);                                               \
    }                                                        \
  }
//---------------------------------------------------------------------------------------
// https://stackoverflow.com/questions/41996828/cuda-reduction-minimum-value-and-index?rq=1
// https://www.olcf.ornl.gov/wp-content/uploads/2019/12/05_Atomics_Reductions_Warp_Shuffle.pdf

__inline__ __device__ void warpReduceMin(int &val, int &idx) {
  for (int offset = warpSize / 2; offset > 0; offset /= 2) {
    int tmpVal = __shfl_down_sync(-1, val, offset);
    int tmpIdx = __shfl_down_sync(-1, idx, offset);
    if (tmpVal < val) {
      val = tmpVal;
      idx = tmpIdx;
    }
  }
}

__inline__ __device__ void blockReduceMin(int &val, int &idx) {
  static __shared__ int values[32];  // Shared mem for 32 partial sums
  static __shared__ int indice[32];  // Shared mem for 32 partial sums

  int lane = threadIdx.x % warpSize;
  int wid = threadIdx.x / warpSize;

  warpReduceMin(val, idx);  // Each warp performs partial reduction

  if (lane == 0) {
    values[wid] = val;  // Write reduced value to shared memory
    indice[wid] = idx;
  }

  __syncthreads();  // Wait for all partial reductions

  // read from shared memory only if that warp existed
  if (threadIdx.x < blockDim.x / warpSize) {
    val = values[lane];
    idx = indice[lane];
  } else {
    val = INT_MAX;
    idx = 0;
  }

  if (wid == 0) warpReduceMin(val, idx);  // Final reduce within first warp
}

__global__ void deviceReduceKernel(int M1, int M2, int *d_cost, int *d_head,
                                   int *d_tail, int *d_PI, int *d_Var,
                                   int mmin) {
  int minVal = mmin;
  int minIdx = -1;

  // set thread ID
  int tid = threadIdx.x;
  int gridSize = blockIdx.x * blockDim.x * 8;
  int idx = M1 + gridSize + tid;

  int e = idx;

  if (idx + 7 * blockDim.x < M2) {
    for (int i = 0; i < 8; i++) {
      e = idx + i * blockDim.x;
      int tmp = d_cost[e] - d_PI[d_head[e]] + d_PI[d_tail[e]];
      if (tmp < minVal) {
        minVal = tmp;
        minIdx = e;
      }
    }
  }

  // synchronize within block
  //__syncthreads();

  //// reduce multiple elements per thread
  // for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < N;
  //     i += blockDim.x * gridDim.x) {
  //  if (in[i] < minVal) {
  //    minVal = in[i];
  //    minIdx = i;  // Added this
  //  }
  //}

  blockReduceMin(minVal, minIdx);

  if (threadIdx.x == 0) {
    // Per la reduction: uso un lock ma lancio il kernel due volte come
    // descritto qui:
    // https://developer.nvidia.com/blog/faster-parallel-reductions-kepler/

    d_Var[blockIdx.x] = minIdx;
  }
}

///-----
__inline__ __device__ void warpReduce(volatile int *vme1, volatile int *vme2,
                                      int tid) {
  if (vme1[tid] > vme1[tid + 32]) {
    vme1[tid] = vme1[tid + 32];
    vme2[tid] = vme2[tid + 32];
  }
  if (vme1[tid] > vme1[tid + 16]) {
    vme1[tid] = vme1[tid + 16];
    vme2[tid] = vme2[tid + 16];
  }
  if (vme1[tid] > vme1[tid + 8]) {
    vme1[tid] = vme1[tid + 8];
    vme2[tid] = vme2[tid + 8];
  }
  if (vme1[tid] > vme1[tid + 4]) {
    vme1[tid] = vme1[tid + 4];
    vme2[tid] = vme2[tid + 4];
  }
  if (vme1[tid] > vme1[tid + 2]) {
    vme1[tid] = vme1[tid + 2];
    vme2[tid] = vme2[tid + 2];
  }
  if (vme1[tid] > vme1[tid + 1]) {
    vme1[tid] = vme1[tid + 1];
    vme2[tid] = vme2[tid + 1];
  }
}

__global__ void fullPricingUnroll88(int M1, int M2, int *d_cost, int *d_head,
                                    int *d_tail, int *d_PI, int *d_Var,
                                    int mmin) {
  __shared__ int viol[BLOCKSIZE];
  __shared__ int best_e[BLOCKSIZE];

  // set thread ID
  int tid = threadIdx.x;
  int gridSize = blockIdx.x * blockDim.x * 8;
  int idx = M1 + gridSize + tid;

  int e = idx;

  viol[tid] = mmin;

  if (idx + 7 * blockDim.x < M2) {
    for (int i = 0; i < 8; i++) {
      e = idx + i * blockDim.x;
      int tmp = d_cost[e] - d_PI[d_head[e]] + d_PI[d_tail[e]];
      if (tmp < viol[tid]) {
        viol[tid] = tmp;
        best_e[tid] = e;
      }
    }
  }

  // synchronize within block
  __syncthreads();

  // in-place reduction in global memory
  if (tid < 512) {
    if (viol[tid] > viol[tid + 512]) {
      viol[tid] = viol[tid + 512];
      best_e[tid] = best_e[tid + 512];
    }
  }
  __syncthreads();
  if (tid < 256) {
    if (viol[tid] > viol[tid + 256]) {
      viol[tid] = viol[tid + 256];
      best_e[tid] = best_e[tid + 256];
    }
  }
  __syncthreads();
  if (tid < 128) {
    if (viol[tid] > viol[tid + 128]) {
      viol[tid] = viol[tid + 128];
      best_e[tid] = best_e[tid + 128];
    }
  }
  __syncthreads();
  if (tid < 64) {
    if (viol[tid] > viol[tid + 64]) {
      viol[tid] = viol[tid + 64];
      best_e[tid] = best_e[tid + 64];
    }
  }
  __syncthreads();

  // unrolling warp
  if (tid < 32) {
    warpReduce(viol, best_e, tid);
    // write result for this block to global mem
    if (tid == 0) {
      int idx = blockIdx.x;
      d_Var[idx] = -1;
      if (viol[tid] < mmin) d_Var[idx] = best_e[tid];
    }
  }
}

#endif