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
   // Base test
   std::string base = "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";
   std::string SEP = "\\";
#ifdef CINECA
   base = "../data/DOTmark/";
   SEP = "/";
#endif

   std::string dtype = "ClassicImages";

   std::string f1 = "data32_1001.csv";
   std::string f2 = "data32_1002.csv";

   Histogram2D h1(base + dtype + SEP + f1);
   Histogram2D h2(base + dtype + SEP + f2);

   fprintf(stdout, "Total weight H1: %ld\n", h1.computeTotalWeight());
   fprintf(stdout, "Total weight H2: %ld\n", h1.computeTotalWeight());

   // Set up configuration option
   Config config;
   config.algo = Algorithm::FlowSimplex;

   config.ground_dist = GroundDistance::L2;
   if (true) { // Exact method
      int n = static_cast<int>(h1.getN());
      // Id using L2, must compute a list of coprimes number
      config.buildCoprimes(n-1);
   }

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

// Default test: Wasserstein distance order 1, ground distance L1, two classic images
void all_dotmark_test() {
   // Base test
   std::string base = "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";
   std::string SEP = "\\";
#ifdef CINECA
   base = "../data/DOTmark/";
   SEP = "/";
#endif

   std::vector<std::string> dirs = { "CauchyDensity", "ClassicImages", "GRFmoderate", "GRFrough", "GRFsmooth",
                                     "LogGRF", "LogitGRF", "MicroscopyImages", "Shapes", "WhiteNoise"
                                   };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = { "32", "64", "128", "256", "512" };

   std::vector<GroundDistance> Gs = { GroundDistance::L1, GroundDistance::Linf, GroundDistance::L2 };

   // Set up configuration option
   Config config;
   config.algo = Algorithm::FlowSimplex;

   for (GroundDistance gd : Gs) {
      config.ground_dist = gd;

      for (const auto& S : Ss) {
         if (config.ground_dist == GroundDistance::L2) { // Exact method
            int n = static_cast<int>(atoi(S.c_str()));
            // Id using L2, must compute a list of coprimes number
            config.buildCoprimes(n - 1);
         }
         for (const auto& dtype : dirs) {
            for (const auto& f11 : Fs) {
               for (const auto& f22 : Fs)
                  if (f11 < f22) {
                     std::string f1 = "data" + S + "_" + f11;
                     std::string f2 = "data" + S + "_" + f22;

                     Histogram2D h1(base + dtype + SEP + f1);
                     Histogram2D h2(base + dtype + SEP + f2);

                     fprintf(stdout, "Total weight H1: %ld\n", h1.computeTotalWeight());
                     fprintf(stdout, "Total weight H2: %ld\n", h1.computeTotalWeight());

                     // Time vars
                     std::chrono::time_point<std::chrono::system_clock> start, end;
                     // Start time.
                     start = std::chrono::system_clock::now();

                     int64_t wd1 = compute_wd1(h1, h2, config);

                     end = std::chrono::system_clock::now();
                     std::chrono::duration<double> inlineTimeElapsed = end - start;
                     std::cout << "RESULT: " << S << ", " << dtype << " " << f1 << " " << f2 << ", "
                               << "Distance: " << wd1
                               << ", Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count() << " ms,"
                               << " Ground: " << config.ground_dist << "\n";

                     start = std::chrono::system_clock::now();
                  }
            }
         }
      }
   }
}

// Main entry point
int main(int argc, char* argv[]) {
   fprintf(stdout, "Running DOTlib v0.2.0\n");

   if (argc == 1) {
      default_test();
      exit(EXIT_SUCCESS);
   }

   if (argc == 2) {
      all_dotmark_test();
      exit(EXIT_SUCCESS);
   }
}
