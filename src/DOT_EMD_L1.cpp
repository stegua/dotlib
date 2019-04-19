/*
*  Main authors:
*     @author Stefano Gualandi <stefano.gualandi@gmail.com>
*
*     @fileoverview Copyright (c) 2017-19, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
*  Last update: April, 2019
*/

#include <chrono>

#include "DOT_KWD_P1.hpp"
#include "emdL1.h"

namespace DOT {

// Compute Kantorovich-Wasserstein Distance (KWD) of order 1, with ground distance L1 (Manhattan)
// using the original code of Link&Okada (see externs/ source files)
void EMD_L1(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg) {

   int m = Mu.size();
   int n = Nu.size();

   size_t s = static_cast<size_t>(sqrt(m));

   std::vector<double> w1;
   w1.reserve(Mu.size());
   for (auto v : Mu.getWs())
      w1.push_back(v);

   std::vector<double> w2;
   w2.reserve(Nu.size());
   for (auto v : Nu.getWs())
      w2.push_back(v);

   // EMD_L1 class (Ling&Okada class implementation)
   EmdL1	em;

   // Timinig output
   auto start = std::chrono::steady_clock::now();

   double fobj = em.EmdDist(&w1[0], &w2[0], s, s);			// 2D EMD_L1

   auto end = std::chrono::steady_clock::now();
   double elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   fprintf(stdout, "%s %d %d Runtime %.6f Value %.f status %d\n",
           msg.c_str(), n, m, elapsed, fobj, 0);

   fflush(stdout);

}

} // End namespace DOT