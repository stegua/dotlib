/**
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#include <chrono>

#include "DOT_BasicTypes.h"
#include "DOT_Histogram2D.h"
#include "DOT_Solvers.h"

using namespace DOT;

// Default test: Wasserstein distance order 1, ground distance L1, two classic images
void default_test() {
   Config config;
   config.ground_dist = GroundDistance::L1;
   config.algo = Algorithm::FlowSimplex;

   // Base test
   std::string base = "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";
   std::string dtype = "ClassicImages";
   std::string f1 = "data32_1001.csv";
   std::string f2 = "data32_1002.csv";

   Histogram2D h1(base + dtype + "\\" + f1);
   Histogram2D h2(base + dtype + "\\" + f2);

   fprintf(stdout, "Total weight H1: %lld\n", h1.computeTotalWeight());
   fprintf(stdout, "Total weight H2: %lld\n", h1.computeTotalWeight());

   // Time vars
   std::chrono::time_point<std::chrono::system_clock> start, end;
   // Start time.
   start = std::chrono::system_clock::now();

   int64_t wd1 = compute_wd1(h1, h2, config);

   end = std::chrono::system_clock::now();
   std::chrono::duration<double> inlineTimeElapsed = end - start;
   std::cout << "RESULT: " << dtype << " " << f1 << " " << f2 << ", "
             << "Distance: " << wd1
             << ", Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count() << " ms \n";

   start = std::chrono::system_clock::now();
}

// Main entry point
int main(int argc, char* argv) {
   fprintf(stdout, "Running DOTlib v0.2.0\n");

   if (argc == 1) {
      default_test();
      exit(EXIT_SUCCESS);
   }
}