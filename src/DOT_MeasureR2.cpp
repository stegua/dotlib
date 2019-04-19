/*
*  Main authors:
*     @author Stefano Gualandi <stefano.gualandi@gmail.com>
*
*     @fileoverview Copyright (c) 2017-19, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
*  Last update: April, 2019
*/

#include <random>
#include <numeric>

#include "DOT_MeasureR2.hpp"

namespace DOT {

// In standarc C++17 it is better to use function std::gcd
int gcd(int _a, int _b) {
   int a = (_a >= 0 ? _a : -_a);
   int b = (_b >= 0 ? _b : -_b);
   while (b != 0) {
      int t = b;
      b = a % b;
      a = t;
   }
   return a;
}

#ifdef _WIN32
const auto& GCD = static_cast<int(*)(int, int)>(std::gcd);
#else
#include <algorithm>
const auto& GCD = gcd;/// static_cast<int(*)(int, int)>(std::__gcd);
#endif

// Create a random measure in R2
MeasureR2 createRandom0N(size_t n, int seed = 13) {
   MeasureR2 mu;
   mu.reserve(n);

   std::random_device rd;  // Will be used to obtain a seed for the random number engine
   std::mt19937 gen(seed); // Standard mersenne_twister_engine seeded with rd()
   std::uniform_real_distribution<> Uniform01(0, n);

   for (size_t i = 0; i < n; i++) {
      PointR2 p = { Uniform01(gen), Uniform01(gen) };
      mu.add(1, p);
   }

   return mu;
}

// Compute co-prime numbers over a lattice n*n with max size L
// List of pair of coprimes number between (-L, L)
std::vector<std::pair<int, int>> DOT::buildCoprimes(int L) {
   std::vector<std::pair<int, int>> coprimes;
   // Reset
   coprimes.clear();
   for (int v = -L; v <= L; ++v)
      for (int w = -L; w <= L; ++w)
         if (!(v == 0 && w == 0) && GCD(v, w) == 1)
            coprimes.emplace_back(v, w);
   // Use as few memory as possible
   coprimes.shrink_to_fit();

   return coprimes;
}

} // End namespace
