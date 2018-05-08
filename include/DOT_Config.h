/**
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "DOT_BasicTypes.h"
#include <numeric>

namespace DOT {

// Type of ground distance
enum GroundDistance { L1=0, L2=1, Linf=2 };

// Convert ground distance enum to string
std::string gd_to_string(GroundDistance gd) {
   if (gd == 0)
      return "L1";
   if (gd == 1)
      return "L2";
   if (gd == 2)
      return "Linf";

   return "Undefined";
}

// Solver algorithm
enum Algorithm { FlowSimplex=0, NetworkSimplexBipartite=1, NetworkSimplexTripartite = 2, CPLEX=3, Gurobi=4 };

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

/**
* @brief Config object
*/
class Config {
 public:
   // Std c'tor
   Config() : seed(13), ground_dist(GroundDistance::L1), algo(Algorithm::FlowSimplex) {}

   // Dump of all parameters
   void dumpParameters() const {
      logger.info("PARAMETERS: %d", seed);
   }

   // Compute the list of coprimes directions as a function of L
   void buildCoprimes(int L) {
      // Reset
      coprimes.clear();
      for (int v = -L; v <= L; ++v)
         for (int w = -L; w <= L; ++w)
            if (!(v == 0 && w == 0) && GCD(v, w) == 1)
               coprimes.emplace_back(v, w);
      // Use as few memory as possible
      coprimes.shrink_to_fit();
   }

   // Seed for random generator device
   int seed;
   // Ground distance
   GroundDistance ground_dist;
   // Type of algorithm for the solver
   Algorithm algo;
   // List of pair of coprimes number between (-L, L)
   std::vector<std::pair<int, int>> coprimes;
};

} // End namespace
