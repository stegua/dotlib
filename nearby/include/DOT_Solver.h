/**
 * @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <omp.h>

#include "DOT_Histogram2D.h"
#include "DOT_NetSimplex.h"
#include "DOT_NetSimplexUnit.h"


struct coprimes_t {
public:
	coprimes_t(int _v, int _w, int _c) : v(_v), w(_w), c_vw(_c) {}
	int v;
	int w;
	int c_vw;
};

typedef std::pair<int, int> int_pair;

struct pair_hash {
	template <class T1, class T2>
	std::size_t operator()(const std::pair<T1, T2>& pair) const {
		return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
	}
};

namespace std {
	template <> struct hash<std::pair<int, int>> {
		inline size_t operator()(const std::pair<int, int>& v) const {
			std::hash<int> int_hasher;
			return int_hasher(v.first) ^ int_hasher(v.second);
		}
	};

} // namespace std

// Zeta block for coordinates vector
#define BLOCKSIZE 1024
#define BLOCKNUM 8

// Taken from: 
#define CHECK(call)                                                            \
{                                                                              \
    const cudaError_t error = call;                                            \
    if (error != cudaSuccess)                                                  \
    {                                                                          \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, "code: %d, reason: %s\n", error,                       \
                cudaGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}


__global__ void naivePricing(
	int Z0,
	int Nv, int* d_V, int* d_W, int* d_Pvw,
	int Npi, int* d_PI,
	int Nvar, int* d_VarB, int* d_VarC)
{
	int tid = threadIdx.x;

	// set thread ID
	__shared__ int viol[BLOCKSIZE];

	viol[tid] = 0;

	if (tid >= Z0)
		return;

	__shared__ int best_c[BLOCKSIZE];
	__shared__ int best_n[BLOCKSIZE];

	int X = blockIdx.x;
	int Y = blockIdx.y;
	int n = gridDim.x;

	int Xv = X + d_V[tid];
	int Yw = Y + d_W[tid];

	//printf("Tid: %d %d %d - Bid: %d %d %d - Bdim: %d %d %d - Gdim: %d %d %d - Pid: %d \n",
	//	threadIdx.x, threadIdx.y, threadIdx.z,
	//	blockIdx.x, blockIdx.y, blockIdx.z,
	//	blockDim.x, blockDim.y, blockDim.z,
	//	gridDim.x, gridDim.y, gridDim.z,
	//	blockDim.x * threadIdx.y + tid);

	if (Xv >= 0 && Xv < n && Yw >= 0 && Yw < n) {
		int p = d_Pvw[tid];
		int idx2 = n * n + (Xv)*n + (Yw);
		viol[tid] = p - d_PI[X * n + Y] + d_PI[idx2];
		best_c[tid] = p;
		best_n[tid] = idx2;
	}

	// synchronize within block
	__syncthreads();

	// in-place reduction in global memory
	//for (int stride = 1; stride < blockDim.x; stride *= 2) {
	//	//	//if ((tid % (2 * stride)) == 0) {
	//	//	//	if (viol[tid] > viol[tid + stride]) {
	//	//	//		viol[tid] = viol[tid + stride];
	//	//	//		best_c[tid] = best_c[tid + stride];
	//	//	//		best_n[tid] = best_n[tid + stride];
	//	//	//	}
	//	//	//}

	//	int index = 2 * stride * tid;

	//	if (index < blockDim.x && viol[index] > viol[index + stride]) {
	//		viol[index] = viol[index + stride];
	//		best_c[index] = best_c[index + stride];
	//		best_n[index] = best_n[index + stride];
	//	}

	//	// synchronize within block
	//	__syncthreads();
	//}

	// in-place reduction in global memory
	for (int stride = blockDim.x / 2; stride > 32; stride >>= 1)
	{
		if (tid < stride && viol[tid] > viol[tid + stride]) {
			viol[tid] = viol[tid + stride];
			best_c[tid] = best_c[tid + stride];
			best_n[tid] = best_n[tid + stride];
		}

		__syncthreads();
	}


	// unrolling warp
	if (tid < 32)
	{
		volatile int* vme1 = viol;
		volatile int* vme2 = best_c;
		volatile int* vme3 = best_n;
		if (vme1[tid] > vme1[tid + 32]) {
			vme1[tid] = vme1[tid + 32];
			vme2[tid] = vme2[tid + 32];
			vme3[tid] = vme3[tid + 32];
		}
		if (vme1[tid] > vme1[tid + 16]) {
			vme1[tid] = vme1[tid + 16];
			vme2[tid] = vme2[tid + 16];
			vme3[tid] = vme3[tid + 16];
		}
		if (vme1[tid] > vme1[tid + 8]) {
			vme1[tid] = vme1[tid + 8];
			vme2[tid] = vme2[tid + 8];
			vme3[tid] = vme3[tid + 8];
		}
		if (vme1[tid] > vme1[tid + 4]) {
			vme1[tid] = vme1[tid + 4];
			vme2[tid] = vme2[tid + 4];
			vme3[tid] = vme3[tid + 4];
		}
		if (vme1[tid] > vme1[tid + 2]) {
			vme1[tid] = vme1[tid + 2];
			vme2[tid] = vme2[tid + 2];
			vme3[tid] = vme3[tid + 2];
		}
		if (vme1[tid] > vme1[tid + 1]) {
			vme1[tid] = vme1[tid + 1];
			vme2[tid] = vme2[tid + 1];
			vme3[tid] = vme3[tid + 1];
		}
		// write result for this block to global mem
		if (tid == 0) {
			int idx = X * n + Y;
			d_VarC[idx] = -1;
			if (viol[tid] < 0) {
				d_VarB[idx] = best_n[tid];
				d_VarC[idx] = best_c[tid];
			}
		}
	}
}


__global__ void naivePricingUnroll2(
	int Nv, int* d_V, int* d_W, int* d_Pvw,
	int* d_PI, int* d_VarB, int* d_VarC)
{
	unsigned int tid = threadIdx.x;

	// set thread ID
	__shared__ int viol[BLOCKSIZE];

	viol[tid] = 0;

	__shared__ int best_c[BLOCKSIZE];
	__shared__ int best_n[BLOCKSIZE];

	int X = blockIdx.x;
	int Y = blockIdx.y;
	int n = gridDim.x;

	int Xv = X + d_V[tid];
	int Yw = Y + d_W[tid];

	//printf("Tid: %d %d %d - Bid: %d %d %d - Bdim: %d %d %d - Gdim: %d %d %d - Pid: %d \n",
	//	threadIdx.x, threadIdx.y, threadIdx.z,
	//	blockIdx.x, blockIdx.y, blockIdx.z,
	//	blockDim.x, blockDim.y, blockDim.z,
	//	gridDim.x, gridDim.y, gridDim.z,
	//	blockDim.x * threadIdx.y + tid);

	if (Xv >= 0 && Xv < n && Yw >= 0 && Yw < n) {
		int p = d_Pvw[tid];
		int idx2 = n * n + (Xv)*n + (Yw);
		viol[tid] = p - d_PI[X * n + Y] + d_PI[idx2];
		best_c[tid] = p;
		best_n[tid] = idx2;
	}

	unsigned int tid2 = blockDim.x * 1 + tid;
	if (tid2 < Nv) {
		int Xv2 = X + d_V[tid2];
		int Yw2 = Y + d_W[tid2];

		if (Xv2 >= 0 && Xv2 < n && Yw2 >= 0 && Yw2 < n) {
			int p2 = d_Pvw[tid2];
			int idx3 = n * n + (Xv2)*n + (Yw2);
			int viol2 = p2 - d_PI[X * n + Y] + d_PI[idx3];
			if (viol2 < viol[tid]) {
				viol[tid] = viol2;
				best_c[tid] = p2;
				best_n[tid] = idx3;
			}
		}
	}

	// synchronize within block
	__syncthreads();

	// in-place reduction in global memory
	//for (int stride = 1; stride < blockDim.x; stride *= 2) {
	//	//if ((tid % (2 * stride)) == 0) {
	//	//	if (viol[tid] > viol[tid + stride]) {
	//	//		viol[tid] = viol[tid + stride];
	//	//		best_c[tid] = best_c[tid + stride];
	//	//		best_n[tid] = best_n[tid + stride];
	//	//	}
	//	//}

	//	int index = 2 * stride * tid;

	//	if (index < blockDim.x && viol[index] > viol[index + stride]) {
	//		viol[index] = viol[index + stride];
	//		best_c[index] = best_c[index + stride];
	//		best_n[index] = best_n[index + stride];
	//	}

	//	// synchronize within block
	//	__syncthreads();
	//}

	// in-place reduction in global memory
	for (int stride = blockDim.x / 2; stride > 32; stride >>= 1)
	{
		if (tid < stride && viol[tid] > viol[tid + stride]) {
			viol[tid] = viol[tid + stride];
			best_c[tid] = best_c[tid + stride];
			best_n[tid] = best_n[tid + stride];
		}

		__syncthreads();
	}

	// unrolling warp
	if (tid < 32)
	{
		volatile int* vme1 = viol;
		volatile int* vme2 = best_c;
		volatile int* vme3 = best_n;
		if (vme1[tid] > vme1[tid + 32]) {
			vme1[tid] = vme1[tid + 32];
			vme2[tid] = vme2[tid + 32];
			vme3[tid] = vme3[tid + 32];
		}
		if (vme1[tid] > vme1[tid + 16]) {
			vme1[tid] = vme1[tid + 16];
			vme2[tid] = vme2[tid + 16];
			vme3[tid] = vme3[tid + 16];
		}
		if (vme1[tid] > vme1[tid + 8]) {
			vme1[tid] = vme1[tid + 8];
			vme2[tid] = vme2[tid + 8];
			vme3[tid] = vme3[tid + 8];
		}
		if (vme1[tid] > vme1[tid + 4]) {
			vme1[tid] = vme1[tid + 4];
			vme2[tid] = vme2[tid + 4];
			vme3[tid] = vme3[tid + 4];
		}
		if (vme1[tid] > vme1[tid + 2]) {
			vme1[tid] = vme1[tid + 2];
			vme2[tid] = vme2[tid + 2];
			vme3[tid] = vme3[tid + 2];
		}
		if (vme1[tid] > vme1[tid + 1]) {
			vme1[tid] = vme1[tid + 1];
			vme2[tid] = vme2[tid + 1];
			vme3[tid] = vme3[tid + 1];
		}
		// write result for this block to global mem
		if (tid == 0) {
			int idx = X * n + Y;
			d_VarC[idx] = -1;
			if (viol[tid] < 0) {
				d_VarB[idx] = best_n[tid];
				d_VarC[idx] = best_c[tid];
			}
		}
	}
}

__global__ void naivePricingUnroll(
	int Nv, int* d_V, int* d_W, int* d_Pvw,
	int* d_PI, int* d_Var)
{

	// set thread ID
	__shared__ int viol[BLOCKSIZE];
	__shared__ int best_c[BLOCKSIZE];
	__shared__ int best_n[BLOCKSIZE];

	int X = blockIdx.x;
	int Y = blockIdx.y;
	int n = gridDim.x;

	//printf("Tid: %d %d %d - Bid: %d %d %d - Bdim: %d %d %d - Gdim: %d %d %d - Pid: %d \n",
	//	threadIdx.x, threadIdx.y, threadIdx.z,
	//	blockIdx.x, blockIdx.y, blockIdx.z,
	//	blockDim.x, blockDim.y, blockDim.z,
	//	gridDim.x, gridDim.y, gridDim.z,
	//	blockDim.x * threadIdx.y + tid);

	unsigned int tid = threadIdx.x;
	unsigned int tid2 = blockDim.x * 1 + tid;
	unsigned int tid3 = blockDim.x * 2 + tid;
	unsigned int tid4 = blockDim.x * 3 + tid;
	unsigned int tid5 = blockDim.x * 4 + tid;
	unsigned int tid6 = blockDim.x * 5 + tid;
	unsigned int tid7 = blockDim.x * 6 + tid;
	unsigned int tid8 = blockDim.x * 7 + tid;

	int Xv = X + d_V[tid];
	int Yw = Y + d_W[tid];

	int Xv2 = X + d_V[tid2];
	int Yw2 = Y + d_W[tid2];

	int Xv3 = X + d_V[tid3];
	int Yw3 = Y + d_W[tid3];

	int Xv4 = X + d_V[tid4];
	int Yw4 = Y + d_W[tid4];

	int Xv5 = X + d_V[tid5];
	int Yw5 = Y + d_W[tid5];

	int Xv6 = X + d_V[tid6];
	int Yw6 = Y + d_W[tid6];

	int Xv7 = X + d_V[tid7];
	int Yw7 = Y + d_W[tid7];

	int Xv8 = X + d_V[tid8];
	int Yw8 = Y + d_W[tid8];

	int NN = n * n;
	int DI = d_PI[X * n + Y];

	viol[tid] = 0;

	if (Xv >= 0 && Xv < n && Yw >= 0 && Yw < n) {
		int p = d_Pvw[tid];
		int idx2 = NN + (Xv)*n + (Yw);
		viol[tid] = p - DI + d_PI[idx2];
		best_c[tid] = p;
		best_n[tid] = idx2;
	}

	if (Xv2 >= 0 && Xv2 < n && Yw2 >= 0 && Yw2 < n) {
		int p2 = d_Pvw[tid2];
		int idx3 = NN + (Xv2)*n + (Yw2);
		int viol2 = p2 - DI + d_PI[idx3];
		if (viol2 < viol[tid]) {
			viol[tid] = viol2;
			best_c[tid] = p2;
			best_n[tid] = idx3;
		}
	}

	if (Xv3 >= 0 && Xv3 < n && Yw3 >= 0 && Yw3 < n) {
		int p3 = d_Pvw[tid3];
		int idx4 = NN + (Xv3)*n + (Yw3);
		int viol3 = p3 - DI + d_PI[idx4];
		if (viol3 < viol[tid]) {
			viol[tid] = viol3;
			best_c[tid] = p3;
			best_n[tid] = idx4;
		}
	}

	if (Xv4 >= 0 && Xv4 < n && Yw4 >= 0 && Yw4 < n) {
		int p4 = d_Pvw[tid4];
		int idx5 = NN + (Xv4)*n + (Yw4);
		int viol4 = p4 - DI + d_PI[idx5];
		if (viol4 < viol[tid]) {
			viol[tid] = viol4;
			best_c[tid] = p4;
			best_n[tid] = idx5;
		}
	}

	if (Xv5 >= 0 && Xv5 < n && Yw5 >= 0 && Yw5 < n) {
		int p = d_Pvw[tid5];
		int idx0 = NN + (Xv5)*n + (Yw5);
		int viol5 = p - DI + d_PI[idx0];
		if (viol5 < viol[tid]) {
			viol[tid] = viol5;
			best_c[tid] = p;
			best_n[tid] = idx0;
		}
	}

	if (Xv6 >= 0 && Xv6 < n && Yw6 >= 0 && Yw6 < n) {
		int p = d_Pvw[tid6];
		int idx0 = NN + (Xv6)*n + (Yw6);
		int viol6 = p - DI + d_PI[idx0];
		if (viol6 < viol[tid]) {
			viol[tid] = viol6;
			best_c[tid] = p;
			best_n[tid] = idx0;
		}
	}

	if (Xv7 >= 0 && Xv7 < n && Yw7 >= 0 && Yw7 < n) {
		int p = d_Pvw[tid7];
		int idx0 = NN + (Xv7)*n + (Yw7);
		int viol7 = p - DI + d_PI[idx0];
		if (viol7 < viol[tid]) {
			viol[tid] = viol7;
			best_c[tid] = p;
			best_n[tid] = idx0;
		}
	}

	if (Xv8 >= 0 && Xv8 < n && Yw8 >= 0 && Yw8 < n) {
		int p = d_Pvw[tid8];
		int idx0 = NN + (Xv8)*n + (Yw8);
		int viol8 = p - DI + d_PI[idx0];
		if (viol8 < viol[tid]) {
			viol[tid] = viol8;
			best_c[tid] = p;
			best_n[tid] = idx0;
		}
	}

	__syncthreads();

	// in-place reduction in global memory
	for (int stride = blockDim.x / 2; stride > 32; stride >>= 1)
	{
		if (tid < stride && viol[tid] > viol[tid + stride]) {
			viol[tid] = viol[tid + stride];
			best_c[tid] = best_c[tid + stride];
			best_n[tid] = best_n[tid + stride];
		}

		__syncthreads();
	}

	// unrolling warp
	if (tid < 32)
	{
		volatile int* vme1 = viol;
		volatile int* vme2 = best_c;
		volatile int* vme3 = best_n;
		if (vme1[tid] > vme1[tid + 32]) {
			vme1[tid] = vme1[tid + 32];
			vme2[tid] = vme2[tid + 32];
			vme3[tid] = vme3[tid + 32];
		}
		if (vme1[tid] > vme1[tid + 16]) {
			vme1[tid] = vme1[tid + 16];
			vme2[tid] = vme2[tid + 16];
			vme3[tid] = vme3[tid + 16];
		}
		if (vme1[tid] > vme1[tid + 8]) {
			vme1[tid] = vme1[tid + 8];
			vme2[tid] = vme2[tid + 8];
			vme3[tid] = vme3[tid + 8];
		}
		if (vme1[tid] > vme1[tid + 4]) {
			vme1[tid] = vme1[tid + 4];
			vme2[tid] = vme2[tid + 4];
			vme3[tid] = vme3[tid + 4];
		}
		if (vme1[tid] > vme1[tid + 2]) {
			vme1[tid] = vme1[tid + 2];
			vme2[tid] = vme2[tid + 2];
			vme3[tid] = vme3[tid + 2];
		}
		if (vme1[tid] > vme1[tid + 1]) {
			vme1[tid] = vme1[tid + 1];
			vme2[tid] = vme2[tid + 1];
			vme3[tid] = vme3[tid + 1];
		}
		// write result for this block to global mem
		if (tid == 0) {
			int idx = X * n + Y;
			d_Var[idx] = -1;
			if (viol[tid] < 0) {
				d_Var[idx] = best_c[tid];
				d_Var[NN + idx] = best_n[tid];
			}
		}
	}
}

namespace DOT {

	// Solver class, which wrapper the Network Simplex algorithm
	class Solver {
	public:
		// Standard c'tor
		Solver()
			: _runtime(0.0), _n_log(0), L(-1), verbosity(DOT_VAL_INFO), recode(""),
			opt_tolerance(0), timelimit(std::numeric_limits<double>::max()) {}

		// Setter/getter for parameters
		std::string getStrParam(const std::string& name) const {
			if (name == DOT_PAR_METHOD)
				return method;
			if (name == DOT_PAR_ALGORITHM)
				return algorithm;
			if (name == DOT_PAR_VERBOSITY)
				return verbosity;
			if (name == DOT_PAR_RECODE)
				return recode;
			return "ERROR getStrParam: wrong parameter ->" + name;
		}

		double getDblParam(const std::string& name) const {
			if (name == DOT_PAR_TIMELIMIT)
				return timelimit;
			if (name == DOT_PAR_OPTTOLERANCE)
				return opt_tolerance;
			return -1;
		}

		void setStrParam(const std::string& name, const std::string& _value) {
			std::string value(_value);
			tolower(value);

			if (name == DOT_PAR_METHOD)
				method = value;

			if (name == DOT_PAR_ALGORITHM)
				algorithm = value;

			if (name == DOT_PAR_VERBOSITY)
				verbosity = value;

			if (name == DOT_PAR_RECODE)
				recode = value;
		}

		void setDblParam(const std::string& name, double value) {
			if (name == DOT_PAR_TIMELIMIT)
				timelimit = value;

			if (name == DOT_PAR_OPTTOLERANCE)
				opt_tolerance = value;
		}

		void dumpParam() const {
			PRINT("Internal parameters: %s %s %s %s %.3f %f %s\n", method.c_str(),
				model.c_str(), algorithm.c_str(), verbosity.c_str(), timelimit,
				opt_tolerance, recode.c_str());
		}

		// Return status of the solver
		std::string status() const {
			if (_status == ProblemType::INFEASIBLE)
				return "Infeasible";
			if (_status == ProblemType::OPTIMAL)
				return "Optimal";
			if (_status == ProblemType::UNBOUNDED)
				return "Unbounded";
			if (_status == ProblemType::TIMELIMIT)
				return "TimeLimit";

			return "Undefined";
		}

		// Return runtime in milliseconds
		double runtime() const { return _runtime; }

		// Number of total iterations of simplex algorithms
		uint64_t iterations() const { return _iterations; }

		// Number of arcs in the model
		uint64_t num_arcs() const { return _num_arcs; }

		// Number of nodes in the model
		uint64_t num_nodes() const { return _num_nodes; }

		//----------------------------------------------------------------------------------------
		// Compute Kantorovich-Wasserstein distance between two measures
		double bipartite(const Histogram2D& A, const Histogram2D& B) {
			int n = A.getN();

			// Compute distances
			std::set<int> tauset;
			for (int v = 0; v < n; ++v)
				for (int w = 0; w < n; ++w)
					tauset.insert(pow(v, 2) + pow(w, 2));

			vector<int> tau;
			for (auto v : tauset)
				tau.push_back(v);

			init_dist_from_to(tau, 0, tau.size() - 1);

			auto start_t = std::chrono::steady_clock::now();

			// Build the graph for min cost flow
			NetSimplex simplex('F', static_cast<int>(2 * n * n),
				static_cast<int>(n * n) * static_cast<int>(n * n));

			// Set the parameters
			simplex.setTimelimit(timelimit);
			simplex.setVerbosity(verbosity);
			simplex.setOptTolerance(opt_tolerance);

			auto ID = [&n](int x, int y) { return x * n + y; };

			// add first d source nodes
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(ID(i, j), A.get(i, j));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(n * n + ID(i, j), -B.get(i, j));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j) {
					for (const auto& p : coprimes) {
						int v = p.v;
						int w = p.w;
						if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
							simplex.addArc(ID(i, j), n * n + ID(i + v, j + w), p.c_vw);
						}
					}
				}

			// Init the simplex
			simplex.run();

			_iterations = simplex.iterations();
			_runtime = simplex.runtime();
			_iterations = simplex.iterations();
			_num_arcs = simplex.num_arcs();
			_num_nodes = simplex.num_nodes();

			auto end_t = std::chrono::steady_clock::now();
			auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
				end_t - start_t)
				.count()) /
				1000;

			double fobj = simplex.totalCost() / A.balance();

			PRINT("BIPARTITE | it: %d, fobj: %.6f, runtime: %.4f (simplex: %.4f), "
				"num_arcs: %ld\n",
				_iterations, fobj, _all, _runtime, _num_arcs);

			return fobj;
		}

		// Compute Kantorovich-Wasserstein distance between two measures
		double tripartite(const Histogram2D& A, const Histogram2D& B) {
			int n = A.getN();

			auto start_t = std::chrono::steady_clock::now();

			// Build the graph for min cost flow
			NetSimplex simplex('F', 3 * n * n, n * n * n);

			// Set the parameters
			simplex.setTimelimit(timelimit);
			simplex.setVerbosity(verbosity);
			simplex.setOptTolerance(opt_tolerance);

			auto ID = [&n](int x, int y) { return x * n + y; };

			// add first d source nodes
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(ID(i, j), A.get(i, j));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(n * n + ID(i, j), 0);

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(2 * n * n + ID(i, j), -B.get(i, j));

			// First layer
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j) {
					for (int h = 0; h < n; ++h) {
						// fprintf(stdout, "(%d, %d)#%d -> (%d, %d)#%d\t", i, j, ID(i, j), h,
						// j,
						//        n * n + ID(h, j));
						simplex.addArc(ID(i, j), n * n + ID(h, j), (int)pow(h - i, 2));
					}
					//        fprintf(stdout, "\n");
				}

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j) {
					for (int h = 0; h < n; ++h) {
						// fprintf(stdout, "(%d, %d)#%d -> (%d, %d)#%d\t", i, j,
						//        n * n + ID(i, j), i, h, 2 * n * n + ID(i, h));
						simplex.addArc(n * n + ID(i, j), 2 * n * n + ID(i, h),
							(int)pow(h - j, 2));
					}
					//       fprintf(stdout, "\n");
				}

			// Init the simplex
			simplex.run();

			_iterations = simplex.iterations();
			_runtime = simplex.runtime();
			_iterations = simplex.iterations();
			_num_arcs = simplex.num_arcs();
			_num_nodes = simplex.num_nodes();

			auto end_t = std::chrono::steady_clock::now();
			auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
				end_t - start_t)
				.count()) /
				1000;

			double fobj = simplex.totalCost() / A.balance();

			PRINT("TRIPARTIE | it: %d, fobj: %.6f, runtime: %.4f (simplex: %.4f), "
				"num_arcs: %ld\n",
				_iterations, fobj, _all, _runtime, _num_arcs);

			return fobj;
		}

		double tripartiteColgen(const Histogram2D& A, const Histogram2D& B) {
			int n = A.getN();
			auto ID = [&n](int x, int y) { return x * n + y; };

			int N = 3 * n * n;
			vector<int> pi(N, 0);

			Vars vars(2 * n * n);
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j) {
					vars[ID(i, j)].a = ID(i, j);
					vars[n * n + ID(i, j)].a = n * n + ID(i, j);
				}

			Vars vnew;
			vnew.reserve(2 * n * n);

			auto start_t = std::chrono::steady_clock::now();

			// Build the graph for min cost flow
			NetSimplex simplex('E', 3 * n * n, 0);

			// Set the parameters
			simplex.setTimelimit(timelimit);
			simplex.setVerbosity(verbosity);
			simplex.setOptTolerance(opt_tolerance);

			// add first d source nodes
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(ID(i, j), A.get(i, j));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(n * n + ID(i, j), 0);

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(2 * n * n + ID(i, j), -B.get(i, j));

			// Init the simplex
			auto _status = simplex.run();

			while (_status != ProblemType::TIMELIMIT) {
				// Take the dual values
				for (int j = 0; j < N; ++j)
					pi[j] = -simplex.potential(j);

				// Solve separation problem:
#pragma omp parallel for collapse(2)
				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j) {
						int best_v1 = 0;
						int best_c1 = -1;
						int best_n1 = 0; // best second node
						int best_v2 = 0;
						int best_c2 = -1;
						int best_n2 = 0; // best second node
						int H1 = ID(i, j);
						int H2 = n * n + ID(i, j);

						for (int h = 0; h < n; ++h) {
							int dist = (h - i) * (h - i);
							int violation = dist - pi[H1] + pi[n * n + ID(h, j)];

							if (violation < best_v1) {
								best_v1 = violation;
								best_c1 = dist;
								best_n1 = n * n + ID(h, j);
							}

							dist = (h - j) * (h - j);
							violation = dist - pi[H2] + pi[2 * n * n + ID(i, h)];

							if (violation < best_v2) {
								best_v2 = violation;
								best_c2 = dist;
								best_n2 = 2 * n * n + ID(i, h);
							}
						}

						vars[H1].b = best_n1;
						vars[H1].c = best_c1;
						vars[H2].b = best_n2;
						vars[H2].c = best_c2;
					}

				// Take all negative reduced cost variables
				vnew.clear();
				for (auto& v : vars) {
					if (v.c > -1)
						vnew.push_back(v);
					v.c = -1;
				}

				if (vnew.empty())
					break;

				std::sort(vnew.begin(), vnew.end(),
					[](const Var& v, const Var& w) { return v.c > w.c; });

				// Replace old constraints with new ones
				int new_arcs = simplex.updateArcs(vnew);

				_status = simplex.reRun();
			}

			_iterations = simplex.iterations();
			_runtime = simplex.runtime();
			_iterations = simplex.iterations();
			_num_arcs = simplex.num_arcs();
			_num_nodes = simplex.num_nodes();

			auto end_t = std::chrono::steady_clock::now();
			auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
				end_t - start_t)
				.count()) /
				1000;

			double fobj = simplex.totalCost() / A.balance();

			PRINT("TRIPARTIE | it: %d, fobj: %.6f, runtime: %.4f (simplex: %.4f), "
				"num_arcs: %ld\n",
				_iterations, fobj, _all, _runtime, _num_arcs);

			return fobj;
		}

		// Compute Kantorovich-Wasserstein distance between two measures
		double phaseTwo(const Histogram2D& A, const Histogram2D& B) {
			int n = A.getN();

			// Compute distances
			std::set<int> tauset;
			for (int v = 0; v < n; ++v)
				for (int w = 0; w < n; ++w)
					tauset.insert(static_cast<int>(pow(v, 2) + pow(w, 2)));

			vector<int> tau;
			for (auto v : tauset)
				tau.push_back(v);
			tau.push_back(tau.back() * 2);

			fprintf(stdout, "distances: %d %d\n", tau.size(), tau[0]);

			int idxL = 0;
			init_dist_from_to(tau, 0, idxL);

			// Build the graph for min cost flow
			NetSimplexUnit simplex('E', static_cast<int>(2 * n * n) + 1, 0);

			auto ID = [&n](int x, int y) { return x * n + y; };
			auto start_t = std::chrono::steady_clock::now();

			{
				// Set the parameters
				simplex.setTimelimit(timelimit);
				simplex.setVerbosity(verbosity);
				simplex.setOptTolerance(opt_tolerance);

				// add first d source nodes
				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j)
						simplex.addNode(ID(i, j), A.get(i, j));

				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j)
						simplex.addNode(n * n + ID(i, j), -B.get(i, j));

				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j) {
						for (const auto& p : coprimes) {
							int v = p.v;
							int w = p.w;
							if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
								simplex.addArc(ID(i, j), n * n + ID(i + v, j + w), p.c_vw);
							}
						}
					}

				int it = 0;
				int64_t fobj = 0;

				// Init the simplex
				simplex.run();
				_iterations = simplex.iterations();

				// Start separation
				while (true) {
					_status = simplex.reRun();

					if (_status == ProblemType::TIMELIMIT)
						break;

					// Check feasibility: if feasible stop
					auto dummyFlow = simplex.dummyFlow();
					fprintf(stdout, "Flow: %ld, idx: %d, tau: %d\n", dummyFlow, idxL,
						tau[idxL]);

					if (abs(dummyFlow) < 1e-06 || idxL >= tau.size())
						break;

					// Add arcs
					init_coprimes(tau[idxL]);
					idxL++;

					for (int i = 0; i < n; ++i)
						for (int j = 0; j < n; ++j) {
							for (const auto& p : coprimes) {
								int v = p.v;
								int w = p.w;
								if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
									simplex.addArc(ID(i, j), n * n + ID(i + v, j + w), p.c_vw);
								}
							}
						}

					++it;
				}

				_runtime = simplex.runtime();
				_iterations = simplex.iterations();
				_num_arcs = simplex.num_arcs();
				_num_nodes = simplex.num_nodes();

				auto end_t = std::chrono::steady_clock::now();
				auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
					end_t - start_t)
					.count()) /
					1000;

				fobj = simplex.totalCost();

				PRINT("it: %d, fobj: %f, all: %f, simplex: %f, num_arcs: %ld\n", it, fobj,
					_all, _runtime, _num_arcs);

				return fobj;
			}

			// ------------------------ PHASE TWO
			// -----------------------------------
			///*idxL = 3;
			// init_coprimes(tau[idxL]);
			// idxL++;*/

			// idxL = 1;
			// init_dist_upto(tau[idxL]);
			// idxL++;

			// fprintf(stdout, "distances: %d %d\n", tau.size(), tau[idxL - 1]);

			//// Build the graph for min cost flow
			// NetSimplex<int, int> simplexTwo('E', static_cast<int>(2 * n * n + 1),
			//	0);

			// for (size_t i = 0; i < n; ++i)
			//	for (size_t j = 0; j < n; ++j)
			//		simplexTwo.addNode(ID(i, j), A.get(i, j));

			// for (size_t i = 0; i < n; ++i)
			//	for (size_t j = 0; j < n; ++j)
			//		simplexTwo.addNode(n * n + ID(i, j), -B.get(i, j));

			// simplexTwo.addNode(2 * n * n, 0);

			// for (size_t i = 0; i < n; ++i)
			//	for (size_t j = 0; j < n; ++j) {
			//		for (const auto& p : coprimes) {
			//			int v = p.v;
			//			int w = p.w;
			//			if (i + v >= 0 && i + v < n && j + w >= 0 && j +
			// w
			//<
			// n) { 				simplexTwo.addArc(ID(i, j), n *
			// n
			// + ID(i
			//+ v,
			// j
			//+ w), p.c_vw);
			//			}
			//		}
			//	}

			NetSimplex simplexTwo(simplex);

			// Arcs to be updated
			vector<size_t> dummy_arcs;
			dummy_arcs.reserve(n * n);
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					dummy_arcs.push_back(simplexTwo.addArc(ID(i, j), 2 * n * n, tau[idxL]));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplexTwo.addArc(2 * n * n, ID(i, j), 0);

			simplexTwo.recomputePotential();

			// Set the parameters
			simplexTwo.setTimelimit(timelimit);
			simplexTwo.setVerbosity(verbosity);
			simplexTwo.setOptTolerance(opt_tolerance);

			_status = simplexTwo.reRun();

			double totF = A.balance();
			while (_status != ProblemType::TIMELIMIT) {
				// Check feasibility: if feasible stop
				auto dummyFlow = simplexTwo.computeDummyFlow(dummy_arcs);
				fprintf(stdout, "Flow: %d, idx: %d, tau: %d, fobj: %.5f\n", dummyFlow,
					idxL, tau[idxL - 1], double(dummyFlow) / totF);

				if (dummyFlow < 1 || idxL >= tau.size() - 1)
					break;

				// Add arcs
				init_coprimes(tau[idxL]);
				idxL++;

				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j) {
						for (const auto& p : coprimes) {
							int v = p.v;
							int w = p.w;
							if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
								simplexTwo.addArc(ID(i, j), n * n + ID(i + v, j + w), p.c_vw);
							}
						}
					}

				simplexTwo.updateArcs(dummy_arcs, tau[idxL]);
				simplexTwo.recomputePotential();

				_status = simplexTwo.reRun();
			}

			_runtime = simplexTwo.runtime();
			_iterations += simplexTwo.iterations();
			_num_arcs = simplexTwo.num_arcs();
			_num_nodes = simplexTwo.num_nodes();

			auto end_t = std::chrono::steady_clock::now();
			auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
				end_t - start_t)
				.count()) /
				1000;

			double fobj = double(simplexTwo.totalCost()) / A.balance();

			PRINT("NEARBY    | it: %d, fobj: %.6f, runtime: %.4f (simplex: %.4f), "
				"num_arcs: %ld\n",
				_iterations, fobj, _all, _runtime, _num_arcs);

			return fobj;
		}

		// Compute Kantorovich-Wasserstein distance between two measures
		double phaseOne(const Histogram2D& A, const Histogram2D& B) {
			int n = A.getN();

			// Compute distances
			std::set<int> tauset;
			for (int v = 0; v < n; ++v)
				for (int w = 0; w < n; ++w)
					tauset.insert(static_cast<int>(pow(v, 2) + pow(w, 2)));

			vector<int> tau;
			for (auto v : tauset)
				tau.push_back(v);

			tau.push_back(tau.back() * 2);

			fprintf(stdout, "distances: %ld %dl\n", tau.size(), tau[0]);

			auto ID = [&n](int x, int y) { return x * n + y; };

			int N = 2 * n * n;
			vector<int> pi(N, 0);

			Vars vars(N);
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					vars[ID(i, j)].a = ID(i, j);

			Vars vnew;
			vnew.reserve(N);

			int TT = 5;
			int idxLO = 0;
			int idxUP = tau.size() / 10;
			init_dist_from_to(tau, idxLO, idxUP);

			// Build the graph for min cost flow
			NetSimplexUnit simplex('E', static_cast<int>(2 * n * n), 0);

			auto start_t = std::chrono::steady_clock::now();

			// Set the parameters
			simplex.setTimelimit(timelimit);
			simplex.setVerbosity(verbosity);
			simplex.setOptTolerance(opt_tolerance);

			// add first d source nodes
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(ID(i, j), A.get(i, j));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(n * n + ID(i, j), -B.get(i, j));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j) {
					for (const auto& p : coprimes) {
						int v = p.v;
						int w = p.w;
						if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
							simplex.addArc(ID(i, j), n * n + ID(i + v, j + w), p.c_vw);
						}
					}
				}

			int64_t fobj = 0;

			// Init the simplex
			simplex.run();

			while (false && _status != ProblemType::TIMELIMIT) {
				// Check feasibility: if feasible stop
				auto dummyFlow = simplex.dummyFlow();
				fprintf(stdout, "Flow: %ld, idx: %d, tau: [%d,%d)\n", dummyFlow, idxUP,
					tau[idxLO], tau[idxUP]);

				if (dummyFlow == 0 || idxUP >= tau.size() - 1)
					break;

				// Add arcs
				idxLO = idxUP;
				idxUP = std::min(idxUP + TT, (int)tau.size());
				init_dist_from_to(tau, idxLO, idxUP);

				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j) {
						for (const auto& p : coprimes) {
							int v = p.v;
							int w = p.w;
							if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
								simplex.addArc(ID(i, j), n * n + ID(i + v, j + w), p.c_vw);
							}
						}
					}
				_status = simplex.reRun();
			}

			_runtime = simplex.runtime();
			_iterations = simplex.iterations();
			_num_arcs = simplex.num_arcs();
			_num_nodes = simplex.num_nodes();

			auto end_t = std::chrono::steady_clock::now();
			auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
				end_t - start_t)
				.count()) /
				1000;

			fobj = simplex.totalCost();

			PRINT("PHASE ONE  | it: %d, fobj: %.6f, runtime: %.4f (simplex: %.4f), "
				"num_arcs: %ld\n",
				_iterations, fobj, _all, _runtime, _num_arcs);

			return fobj;
		}

		double colgenOld(const Histogram2D& A, const Histogram2D& B, int idxL,
			const std::string& msg) {
			int n = A.getN();

			// Compute distances
			std::set<int> tauset;
			for (int v = 0; v < n; ++v)
				for (int w = 0; w < n; ++w)
					tauset.insert(static_cast<int>(pow(v, 2) + pow(w, 2)));

			vector<int> tau;
			for (auto v : tauset)
				tau.push_back(v);
			fprintf(stdout, "distances: %d\n", tau.size());

			//    tau.push_back(tau.back() * 2);

			//    size_t idxL =;
			int TT = std::min<int>(1024, static_cast<int>(tau.size() - 1));
			init_dist_from_to(tau, 0, TT);

			fprintf(stdout, "coprimes size: %d\n", coprimes.size());

			//coprimes.resize(128);

			auto ID = [&n](int x, int y) { return x * n + y; };

			int N = 2 * n * n;
			vector<int> pi(N + 1, 0);

			Vars vars(N + 1);
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					vars[ID(i, j)].a = ID(i, j);
			vars[N].a = N;

			Vars vnew;
			vnew.reserve(N + 1);

			auto start_t = std::chrono::steady_clock::now();

			// Build the graph for min cost flow
			NetSimplex simplex('E', N + 1, 0);

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(ID(i, j), A.get(i, j));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(n * n + ID(i, j), -B.get(i, j));

			simplex.addNode(N, 0);

			vector<size_t> left_arcs, right_arcs;
			left_arcs.reserve(n * n);
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					left_arcs.push_back(simplex.addArc(ID(i, j), N, tau[TT]));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					right_arcs.push_back(simplex.addArc(N, n * n + ID(i, j), 0));

			// Set the parameters
			simplex.setTimelimit(timelimit);
			simplex.setVerbosity(verbosity);
			simplex.setOptTolerance(opt_tolerance);

			_status = simplex.run();

			while (_status != ProblemType::TIMELIMIT) {
				// Take the dual values
				for (int j = 0; j < N; ++j)
					pi[j] = -simplex.potential(j);

				// Solve separation problem:
//#pragma omp parallel for collapse(2)
				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j) {
						int best_v = 0;
						int best_c = -1;
						int best_n = 0; // best second node
						int h = ID(i, j);
						//for (const auto& p : coprimes) 
						for (int PP = 0; PP < std::min<int>(coprimes.size(), BLOCKSIZE * BLOCKNUM); PP++)
						{
							const auto& p = coprimes[PP];
							int v = p.v;
							int w = p.w;
							if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
								int violation = p.c_vw - pi[h] + pi[n * n + ID(i + v, j + w)];
								if (violation < best_v) {
									best_v = violation;
									best_c = p.c_vw;
									best_n = n * n + ID(i + v, j + w);
								}
							}
						}

						// Store most violated cuts for element i
						vars[h].b = best_n;
						vars[h].c = best_c;
					}

				// Take all negative reduced cost variables
				vnew.clear();
				for (auto& v : vars) {
					fprintf(stdout, "%d %d\t", v.b, v.c);
					if (v.c > -1)
						vnew.push_back(v);
					v.c = -1;
				}
				fprintf(stdout, "\n");
				fflush(stdout);

				if (vnew.empty())
					break;

				std::sort(vnew.begin(), vnew.end(),
					[](const Var& v, const Var& w) { return v.c > w.c; });

				// Replace old constraints with new ones
				int new_arcs = simplex.addArcs(vnew);

				_status = simplex.reRun();
			}

			_runtime = simplex.runtime();
			_iterations = simplex.iterations();
			_num_arcs = simplex.num_arcs();
			_num_nodes = simplex.num_nodes();

			auto end_t = std::chrono::steady_clock::now();
			auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
				end_t - start_t)
				.count()) /
				1000;

			double fobj = double(simplex.totalCost()) / A.balance();

			// Upper bound on the missed mass
			auto unmoved = 0.0;// double(simplex.computeDummyFlow(left_arcs)) / A.balance();
			double delta = 0.0;

			vector<double> Aflow(n * n, 0.0);
			vector<double> Bflow(n * n, 0.0);

			Histogram2D AA(A);
			Histogram2D BB(B);

			if (unmoved > 0) {
				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j) {
						Aflow[ID(i, j)] = fabs(simplex.arcFlow(left_arcs[ID(i, j)]));
						AA.set(i, j, Aflow[ID(i, j)]);
						Aflow[ID(i, j)] = Aflow[ID(i, j)] / A.balance();

						Bflow[ID(i, j)] = fabs(simplex.arcFlow(right_arcs[ID(i, j)]));
						BB.set(i, j, Bflow[ID(i, j)]);
						Bflow[ID(i, j)] = Bflow[ID(i, j)] / A.balance();
					}

				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j)
						for (int v = 0; v < n; ++v)
							for (int w = 0; w < n; ++w)
								delta +=
								std::max(0, static_cast<int>(pow(i - v, 2) + pow(j - w, 2)) -
									tau[TT + 1]) *
								(Aflow[ID(i, j)]) * (Bflow[ID(v, w)]) / unmoved;
			}

			// delta = findUB(AA, BB) / A.balance();

			PRINT("COLGEN %s it %lld LB %.6f UB %.6f runtime %.4f simplex %.4f "
				"num_arcs %lld idx %d tau %d maxtau %d residual %.6f\n",
				msg.c_str(), _iterations, fobj, delta, _all, _runtime, _num_arcs, TT,
				tau[TT], tau[tau.size() - 1], unmoved);

			return fobj;
		}

		double colgen(const Histogram2D& A, const Histogram2D& B, int idxL,
			const std::string& msg) {
			int n = A.getN();

			auto ID = [&n](int x, int y) { return x * n + y; };

			int N = 2 * n * n;
			vector<int> pi(N, 0);

			Vars vars(n * n);
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					vars[ID(i, j)].a = ID(i, j);

			Vars vnew;
			vnew.reserve(n * n);

			// Build the graph for min cost flow
			NetSimplex simplex('E', N, 0);

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(ID(i, j), A.get(i, j));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(n * n + ID(i, j), -B.get(i, j));

			// Set the parameters
			simplex.setTimelimit(timelimit);
			simplex.setVerbosity(verbosity);
			simplex.setOptTolerance(opt_tolerance);

			_status = simplex.run();

			auto start_t = std::chrono::steady_clock::now();


			while (_status != ProblemType::TIMELIMIT) {
				// Take the dual values
				for (int j = 0; j < N; ++j)
					pi[j] = -simplex.potential(j);

				// Solve separation problem:
#pragma omp parallel for collapse(2)
				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j) {
						int best_v = 0;
						int best_c = -1;
						int best_n = 0; // best second node
						int h = ID(i, j);
						//for (const auto& p : coprimes) 
						for (int PP = 0; PP < std::min<int>(coprimes.size(), BLOCKNUM * BLOCKSIZE); PP++)
						{
							const auto& p = coprimes[PP];
							int v = p.v;
							int w = p.w;
							if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
								int violation = p.c_vw - pi[h] + pi[n * n + ID(i + v, j + w)];
								if (violation < best_v) {
									best_v = violation;
									best_c = p.c_vw;
									best_n = n * n + ID(i + v, j + w);
								}
							}
						}

						// Store most violated cuts for element i
						vars[h].b = best_n;
						vars[h].c = best_c;
					}

				// Take all negative reduced cost variables
				vnew.clear();
				for (auto& v : vars) {
					if (v.c > -1)
						vnew.push_back(v);
					v.c = -1;
				}

				if (vnew.empty())
					break;

				std::sort(vnew.begin(), vnew.end(),
					[](const Var& v, const Var& w) { return v.c > w.c; });

				// Replace old constraints with new ones
				int new_arcs = simplex.updateArcs(vnew);

				_status = simplex.reRun();
			}

			_runtime = simplex.runtime();
			_iterations = simplex.iterations();
			_num_arcs = simplex.num_arcs();
			_num_nodes = simplex.num_nodes();

			auto end_t = std::chrono::steady_clock::now();
			auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
				end_t - start_t)
				.count()) /
				1000;

			double fobj = double(simplex.totalCost()) / A.balance();

			// Upper bound on the missed mass
			PRINT("COLGEN %s it %lld LB %.6f runtime %.4f simplex %.4f "
				"num_arcs %lld\n",
				msg.c_str(), _iterations, fobj, _all, _runtime, _num_arcs);

			return fobj;
		}

		//------------------------------------------------------------------------------------------
		double colgenCuda(const Histogram2D& A, const Histogram2D& B, int idxL,
			const std::string& msg) {
			int n = A.getN();
			int NN = n * n;
			int NN2 = 2 * n * n;

			auto ID = [&n](int x, int y) { return x * n + y; };

			vector<int> vars(NN, 0);
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					vars[ID(i, j)] = ID(i, j);

			Vars vnew;
			vnew.reserve(NN);

			int* h_V = (int*)malloc(BLOCKNUM * BLOCKSIZE * sizeof(int));
			int* h_W = (int*)malloc(BLOCKNUM * BLOCKSIZE * sizeof(int));
			int* h_Pvw = (int*)malloc(BLOCKNUM * BLOCKSIZE * sizeof(int));

			for (int hh = 0; hh < BLOCKNUM * BLOCKSIZE; hh++)
				if (hh < static_cast<int>(coprimes.size())) {
					auto cc = coprimes[hh];
					h_V[hh] = cc.v;
					h_W[hh] = cc.w;
					h_Pvw[hh] = cc.c_vw;
				}
				else {
					h_V[hh] = 2 * n;
					h_W[hh] = 2 * n;
					h_Pvw[hh] = INT_MAX;
				}

			int Z0 = std::min<int>(BLOCKSIZE * BLOCKNUM, static_cast<int>(coprimes.size()));

			// Data for CUDA support 
			// set up device
			int dev = 0;
			cudaDeviceProp deviceProp;
			cudaGetDeviceProperties(&deviceProp, dev);
			cudaSetDevice(dev);

			int* d_V;
			int* d_W;
			int* d_Pvw;
			size_t Nv = BLOCKNUM * BLOCKSIZE * sizeof(int);
			cudaMalloc((void**)&d_V, Nv);
			cudaMalloc((void**)&d_W, Nv);
			cudaMalloc((void**)&d_Pvw, Nv);
			// Copy once for all
			CHECK(cudaMemcpy(d_V, h_V, Nv, cudaMemcpyHostToDevice));
			CHECK(cudaMemcpy(d_W, h_W, Nv, cudaMemcpyHostToDevice));
			CHECK(cudaMemcpy(d_Pvw, h_Pvw, Nv, cudaMemcpyHostToDevice));

			// Size for dual variables
			size_t Npi = NN2 * sizeof(int);
			int* h_PI;
			int* d_PI;
			CHECK(cudaMallocHost((void**)&h_PI, Npi));
			CHECK(cudaMalloc((void**)&d_PI, Npi));


			size_t Nvar = NN2 * sizeof(int);
			int* h_Var;
			int* d_Var;
			CHECK(cudaMallocHost((void**)&h_Var, Nvar));
			CHECK(cudaMalloc((void**)&d_Var, Nvar));

			const dim3 blockSize(BLOCKSIZE, 1);
			const dim3 gridSize(n, n, 1);


			// Build the graph for min cost flow
			auto start_t = std::chrono::steady_clock::now();

			NetSimplex simplex('E', NN2, 0);

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(ID(i, j), A.get(i, j));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(n * n + ID(i, j), -B.get(i, j));

			// Set the parameters
			simplex.setTimelimit(timelimit);
			simplex.setVerbosity(verbosity);
			simplex.setOptTolerance(opt_tolerance);

			_status = simplex.run();

			while (_status != ProblemType::TIMELIMIT) {
				// Take the dual values
				for (int j = 0; j < NN2; ++j)
					h_PI[j] = -simplex.potential(j);

				CHECK(cudaMemcpy(d_PI, h_PI, Npi, cudaMemcpyHostToDevice));

				naivePricingUnroll << <gridSize, blockSize >> > (Z0, d_V, d_W, d_Pvw, d_PI, d_Var);

				CHECK(cudaMemcpy(h_Var, d_Var, Nvar, cudaMemcpyDeviceToHost));

				// Take all negative reduced cost variables
				vnew.clear();
				for (int h = 0, h_max = NN; h < h_max; ++h)
					if (h_Var[h] > -1)
						vnew.emplace_back(vars[h], h_Var[NN + h], h_Var[h]);

				if (vnew.empty())
					break;

				std::sort(vnew.begin(), vnew.end(),
					[](const Var& v, const Var& w) { return v.c > w.c; });

				// Replace old constraints with new ones
				int new_arcs = simplex.updateArcs(vnew);

				_status = simplex.reRun();
			}

			_runtime = simplex.runtime();
			_iterations = simplex.iterations();
			_num_arcs = simplex.num_arcs();
			_num_nodes = simplex.num_nodes();

			auto end_t = std::chrono::steady_clock::now();
			auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
				end_t - start_t)
				.count()) /
				1000;

			double fobj = double(simplex.totalCost()) / A.balance();


			// Upper bound on the missed mass
			double delta = 0.0;

			PRINT("COCUDA %s it %lld LB %.6f UB %.6f runtime %.4f simplex %.4f "
				"num_arcs %lld\n",
				msg.c_str(), _iterations, fobj, delta, _all, _runtime, _num_arcs);

			// free device memory
			cudaFree(d_V);
			cudaFree(d_W);
			cudaFree(d_Pvw);
			cudaFree(d_PI);
			cudaFree(d_Var);

			cudaFreeHost(h_Var);
			cudaFreeHost(h_PI);

			free(h_V);
			free(h_W);
			free(h_Pvw);
			cudaDeviceReset();

			return fobj;
		}

		double nearbyUB(const Histogram2D& A, const Histogram2D& B, int idxL,
			const std::string& msg) {
			int n = A.getN();

			auto ID = [&n](int x, int y) { return x * n + y; };

			int N = 2 * n * n;
			vector<int> pi(N, 0);

			Vars vars(N);
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					vars[ID(i, j)].a = ID(i, j);

			Vars vnew;
			vnew.reserve(N);

			auto start_t = std::chrono::steady_clock::now();

			// Build the graph for min cost flow
			NetSimplex simplex('E', N, 0);

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(ID(i, j), A.get(i, j));

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					simplex.addNode(n * n + ID(i, j), -B.get(i, j));

			// Set the parameters
			simplex.setTimelimit(timelimit);
			simplex.setVerbosity(verbosity);
			simplex.setOptTolerance(opt_tolerance);

			_status = simplex.run();

			while (_status != ProblemType::TIMELIMIT) {
				// Take the dual values
				for (int j = 0; j < N; ++j)
					pi[j] = -simplex.potential(j);

				// Solve separation problem:
#pragma omp parallel for collapse(2)
				for (int i = 0; i < n; ++i)
					for (int j = 0; j < n; ++j) {
						int best_v = 0;
						int best_c = -1;
						int best_n = 0; // best second node
						int h = ID(i, j);
						for (const auto& p : coprimes) {
							int v = p.v;
							int w = p.w;
							if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
								int violation = p.c_vw - pi[h] + pi[n * n + ID(i + v, j + w)];
								if (violation < best_v) {
									best_v = violation;
									best_c = p.c_vw;
									best_n = n * n + ID(i + v, j + w);
								}
							}
						}

						// Store most violated cuts for element i
						vars[h].b = best_n;
						vars[h].c = best_c;
					}

				// Take all negative reduced cost variables
				vnew.clear();
				for (auto& v : vars) {
					if (v.c > -1)
						vnew.push_back(v);
					v.c = -1;
				}

				if (vnew.empty())
					break;

				std::sort(vnew.begin(), vnew.end(),
					[](const Var& v, const Var& w) { return v.c > w.c; });

				// Replace old constraints with new ones
				int new_arcs = simplex.updateArcs(vnew);

				_status = simplex.reRun();
			}

			_runtime = simplex.runtime();
			_iterations = simplex.iterations();
			_num_arcs = simplex.num_arcs();
			_num_nodes = simplex.num_nodes();

			auto end_t = std::chrono::steady_clock::now();
			auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
				end_t - start_t)
				.count()) /
				1000;

			double fobj = double(simplex.totalCost()) / A.balance();

			// Upper bound on the missed mass
			// double(simplex.computeDummyFlow(left_arcs)) / A.balance();
			double delta = 0.0;

			PRINT("COLEGN %s it %lld LB %.6f UB %.6f runtime %.4f simplex %.4f "
				"num_arcs %lld\n",
				msg.c_str(), _iterations, fobj, delta, _all, _runtime, _num_arcs);

			return fobj;
		}

#ifdef __MY_LEMON__
		double Winf(const Histogram2D& A, const Histogram2D& B) {
			int s = A.getN();
			int d = s * s;

			// Compute distances
			std::set<int> tauset;
			for (int v = 0; v < s; ++v)
				for (int w = 0; w < s; ++w)
					tauset.insert(pow(v, 2) + pow(w, 2));

			vector<int> tau;
			for (auto v : tauset)
				tau.push_back(v);

			auto start_t = std::chrono::steady_clock::now();

			int idxL = tau.size() / 10;
			init_dist_from_to(tau, tau[0], tau[idxL]);
			idxL++;

			fprintf(stdout, "distances: %d %d\n", tau.size(), tau[idxL - 1]);

			LemonGraph g;

			auto ID = [&s](int x, int y) { return x * s + y; };

			// add d nodes for each histrogam (d+1) source, (d+2) target
			std::vector<LemonGraph::Node> nodes;
			// add first d source nodes
			for (int i = 0; i < s; ++i)
				for (int j = 0; j < s; ++j)
					nodes.emplace_back(g.addNode());

			for (int i = 0; i < s; ++i)
				for (int j = 0; j < s; ++j)
					nodes.emplace_back(g.addNode());

			// Sink node
			nodes.emplace_back(g.addNode());
			nodes.emplace_back(g.addNode());

			auto S = nodes[nodes.size() - 2];
			auto T = nodes[nodes.size() - 1];

			std::vector<LemonGraph::Arc> arcs;
			std::vector<int> a_cap;

			for (int i = 0; i < s; ++i)
				for (int j = 0; j < s; ++j) {
					for (const auto& p : coprimes) {
						int v = p.v;
						int w = p.w;
						if (i + v >= 0 && i + v < s && j + w >= 0 && j + w < s) {
							arcs.emplace_back(
								g.addArc(nodes[ID(i, j)], nodes[s * s + ID(i + v, j + w)]));
							a_cap.emplace_back(std::min(A.get(i, j), B.get(i + v, j + w)));
						}
					}
				}

			for (int i = 0; i < s; ++i)
				for (int j = 0; j < s; ++j) {
					arcs.emplace_back(g.addArc(S, nodes[ID(i, j)]));
					a_cap.emplace_back(A.get(i, j));
				}

			for (int i = 0; i < s; ++i)
				for (int j = 0; j < s; ++j) {
					arcs.emplace_back(g.addArc(nodes[s * s + ID(i, j)], T));
					a_cap.emplace_back(B.get(i, j));
				}

			fprintf(stdout, "Input graph created with %d nodes and %d arcs\n",
				countNodes(g), countArcs(g));
			ListDigraph::ArcMap<int> u_i(g);
			for (int i = 0, i_max = arcs.size(); i < i_max; ++i) {
				const auto& a = arcs[i];
				u_i[a] = a_cap[i];
			}

			Preflow<LemonGraph> solver(g, u_i, S, T);
			solver.runMinCut();

			auto end_t = std::chrono::steady_clock::now();
			auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
				end_t - start_t)
				.count()) /
				1000;

			fprintf(stdout, "max flow value: %d, time: %.4f\n", solver.flowValue(),
				_all);
		}

		// Compute Kantorovich-Wasserstein distance between two measures
		double lemon(const Histogram2D& A, const Histogram2D& B) {
			int s = A.getN();
			int d = s * s;

			// Compute distances
			std::set<int> tauset;
			for (int v = 0; v < s; ++v)
				for (int w = 0; w < s; ++w)
					tauset.insert(pow(v, 2) + pow(w, 2));

			vector<int> tau;
			for (auto v : tauset)
				tau.push_back(v);

			auto start_t = std::chrono::steady_clock::now();

			int idxL = 3;
			init_dist_upto(tau[idxL]);
			idxL++;

			fprintf(stdout, "distances: %d %d\n", tau.size(), tau[idxL - 1]);

			LemonGraph g;

			auto ID = [&s](int x, int y) { return x * s + y; };

			// add d nodes for each histrogam (d+1) source, (d+2) target
			std::vector<LemonGraph::Node> nodes;
			// add first d source nodes
			for (int i = 0; i < s; ++i)
				for (int j = 0; j < s; ++j)
					nodes.emplace_back(g.addNode());

			for (int i = 0; i < s; ++i)
				for (int j = 0; j < s; ++j)
					nodes.emplace_back(g.addNode());

			// Sink node
			nodes.emplace_back(g.addNode());

			std::vector<LemonGraph::Arc> arcs;
			std::vector<int64_t> a_costs;
			std::vector<int64_t> a_cap;

			for (int i = 0; i < s; ++i)
				for (int j = 0; j < s; ++j) {
					for (const auto& p : coprimes) {
						int v = p.v;
						int w = p.w;
						if (i + v >= 0 && i + v < s && j + w >= 0 && j + w < s) {
							arcs.emplace_back(
								g.addArc(nodes[ID(i, j)], nodes[s * s + ID(i + v, j + w)]));
							a_costs.emplace_back(-tau[idxL] + p.c_vw);
							a_cap.emplace_back(std::min(A.get(i, j), B.get(i + v, j + w)));
						}
					}
				}

			for (int i = 0; i < s; ++i)
				for (int j = 0; j < s; ++j) {
					arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[2 * s * s]));
					a_costs.emplace_back(0);
					a_cap.emplace_back(A.get(i, j));
				}

			for (int i = 0; i < s; ++i)
				for (int j = 0; j < s; ++j) {
					arcs.emplace_back(g.addArc(nodes[2 * s * s], nodes[s * s + ID(i, j)]));
					a_costs.emplace_back(0);
					a_cap.emplace_back(B.get(i, j));
				}

			fprintf(stdout, "Input graph created with %d nodes and %d arcs\n",
				countNodes(g), countArcs(g));

			LemonSimplex simplex(g);

			// lower and upper bounds, cost
			ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
			ListDigraph::ArcMap<LimitValueType> c_i(g);

			// FLow balance
			ListDigraph::NodeMap<LimitValueType> b_i(g);
			{
				int idx = 0;
				for (int i = 0; i < s; ++i)
					for (int j = 0; j < s; ++j)
						b_i[nodes[idx++]] = LimitValueType(A.get(i, j));

				for (int i = 0; i < s; ++i)
					for (int j = 0; j < s; ++j)
						b_i[nodes[idx++]] = LimitValueType(-B.get(i, j));

				b_i[nodes[idx++]] = LimitValueType(0);
			}

			// Add all edges
			for (int i = 0, i_max = arcs.size(); i < i_max; ++i) {
				const auto& a = arcs[i];
				l_i[a] = 0;
				u_i[a] = a_cap[i];
				c_i[a] = a_costs[i];
			}

			// set lower/upper bounds, cost
			simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

			// simplex.supplyType(NetworkSimplex<LemonGraph, LimitValueType,
			// LimitValueType>::LEQ);

			// Solve the problem to compute the distance
			NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType
				ret = simplex.run();

			switch (ret) {
			case NetworkSimplex<LemonGraph>::INFEASIBLE:
				fprintf(stdout, "NetworkSimplex<Graph>::INFEASIBLE\n");
				break;
			case NetworkSimplex<LemonGraph>::OPTIMAL:
				fprintf(stdout, "NetworkSimplex<Graph>::OPTIMAL\n");
				break;
			case NetworkSimplex<LemonGraph>::UNBOUNDED:
				fprintf(stdout, "NetworkSimplex<Graph>::UNBOUNDED\n");
				break;
			}

			double ff = 0;
			for (int tt = arcs.size() - s * s, tt_max = arcs.size(); tt < tt_max;
				++tt) {
				ff += simplex.flow(arcs[tt]);
				/*if (simplex.flow(arcs[tt]) > 0)*/
				fprintf(stdout, "%d\t", simplex.flow(arcs[tt]));
			}

			fprintf(stdout, "transported %.f %.f\n", ff, (double)A.balance());

			double fobj = tau[idxL] + double(simplex.totalCost()) / A.balance();

			PRINT("LEMON    | it: %d, fobj: %.6f, tau: %d, simplex: %.4f, num_arcs: "
				"%ld\n",
				_iterations, fobj, tau[idxL - 1], _runtime, arcs.size());

			return fobj;
		}
#endif

		// Upper bound: find a feasible plane
		double findUB(const Histogram2D& _A, const Histogram2D& _B) {
			Histogram2D A(_A);
			Histogram2D B(_B);

			int n = A.getN();

			// pairs ordered
			init_all(n, true);

			double delta = 0;

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j) {
					int Aij = A.get(i, j);
					for (const auto& p : coprimes) {
						int v = p.v;
						int w = p.w;
						if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
							int Bij = B.get(i + v, j + w);
							if (Aij <= Bij) {
								delta += (double)p.c_vw * Aij;
								B.subtract(i + v, j + w, Aij);
								A.subtract(i, j, Aij);
								break;
							}
							else {
								delta += (double)p.c_vw * Bij;
								Aij = Aij - Bij;
								B.subtract(i + v, j + w, Bij);
								A.subtract(i, j, Bij);
							}
						}
					}
				}

			// fprintf(stdout, "balance A: %lld, balance B: %lld, value: %.6f\n",
			//        A.balance(), B.balance(), delta / _A.balance());
			return delta;
		}

		//--------------------------------------------------------------------------
		void init_coprimes(int L) {
			coprimes.clear();
			for (int v = -L; v <= L; ++v)
				for (int w = -L; w <= L; ++w)
					if (pow(v, 2) + pow(w, 2) == L)
						coprimes.emplace_back(v, w, L);
			coprimes.shrink_to_fit();
		}

		//--------------------------------------------------------------------------
		void init_dist_from_to(const vector<int>& data, int LO, int UP,
			bool order = false) {
			coprimes.clear();
			for (int h = LO; h < UP; h++) {
				int L = data[h];
				for (int v = -L; v <= L; ++v)
					for (int w = -L; w <= L; ++w)
						if (pow(v, 2) + pow(w, 2) == L)
							coprimes.emplace_back(v, w, L);
			}
			if (order)
				std::sort(coprimes.begin(), coprimes.end(),
					[](const auto& a, const auto& b) { return a.c_vw < b.c_vw; });

			//fprintf(stdout, "__device__ int s_V[8192] = {");
			//for (int i = 0; i < 1024 * 8; i++)
			//	fprintf(stdout, "%d, ", coprimes[i].v);
			//fprintf(stdout, "};\n");
			//fprintf(stdout, "__device__ int s_W[8192] = {");
			//for (int i = 0; i < 1024 * 8; i++)
			//	fprintf(stdout, "%d, ", coprimes[i].w);
			//fprintf(stdout, "};\n");
			//fprintf(stdout, "__device__ int s_Pvw[8192] = {");
			//for (int i = 0; i < 1024 * 8; i++)
			//	fprintf(stdout, "%d, ", coprimes[i].c_vw);
			//fprintf(stdout, "};\n");

			coprimes.shrink_to_fit();
		}

		//--------------------------------------------------------------------------
		void init_all(int n, bool order = false) {
			coprimes.clear();
			for (int v = -n + 1; v < n; ++v)
				for (int w = -n + 1; w < n; ++w)
					coprimes.emplace_back(v, w, pow(v, 2) + pow(w, 2));
			coprimes.shrink_to_fit();

			if (order)
				std::sort(coprimes.begin(), coprimes.end(),
					[](const auto& a, const auto& b) { return a.c_vw < b.c_vw; });
		}

		// List of pair of coprimes number between (-L, L)
		std::vector<coprimes_t> coprimes;

	private:
		// Status of the solver
		ProblemType _status;

		// Runtime in milliseconds
		double _runtime;

		// Number of iterations
		uint64_t _iterations;
		uint64_t _num_nodes;
		uint64_t _num_arcs;

		// Interval for logging iterations in the simplex algorithm
		// (if _n_log=0 no logs at all)
		int _n_log;

		// Approximation parameter
		int L;


		// Method to solve the problem
		std::string method;
		// Model to solve the problem
		std::string model;
		// Algorithm to solve the corresponding problem
		std::string algorithm;
		// Verbosity of the log
		std::string verbosity;
		// Recode the coordinates as consecutive integers
		std::string recode;
		// Tolerance for pricing
		double opt_tolerance;
		// Time limit for runtime of the algorithm
		double timelimit;

	}; // namespace DOT

} // namespace DOT
