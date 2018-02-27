/**
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "DOT_Histogram2D.h"
#include "DOT_Config.h"
#include "DOT_WD1.h"

namespace DOT {

/**
* @brief Compute Wasserstein distance of order p
*/
int64_t compute_wdp(const Histogram2D& h1, const Histogram2D& h2, const Config& config) {

   return 0;
}

/**
* @brief Compute Wasserstein distance of order 1
*/
int64_t compute_wd1(const Histogram2D& h1, const Histogram2D& h2, const Config& config) {
   if (h1.getN() != h2.getN()) {
      logger.error("ERROR: The two histograms must have the same dimension!");
      return -1;
   }

   if (config.algo == Algorithm::FlowSimplex) {
      if (config.ground_dist == GroundDistance::L1)
         return solve_network_L1(h1, h2);
      if (config.ground_dist == GroundDistance::Linf)
         return solve_network_Linf(h1, h2);
      if (config.ground_dist == GroundDistance::L2) {
         return solve_network_L2(h1, h2, config.coprimes);
      }
   }
   return -1;
}

};
