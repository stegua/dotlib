/**
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "DOT_BasicTypes.h"

namespace DOT {

// Type of ground distance
enum GroundDistance { L1=0, L2=1, Linf=2 };

// Solver algorithm
enum Algorithm { FlowSimplex=0, CPLEX=1, Gurobi=2 };


/**
* @brief Config object
*/
class Config {
 public:
   // Std c'tor
   Config() : seed(13), ground_dist(GroundDistance::L1), algo(Algorithm::FlowSimplex) {}

   // Dump of all parameters
   void dumpParameters() const {
      fprintf(stdout, "PARAMETERS: %d\n", seed);
   }

   // Seed for random generator device
   int seed;

   // Ground distance
   GroundDistance ground_dist;

   // Type of algorithm for the solver
   Algorithm algo;
};

} // End namespace