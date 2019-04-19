/*
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#include <chrono>

// Args parser library
#include <args.hxx>

// Yocta Logger library
#include <yocta_logger.hh>

// Global variable logger
namespace DOT {
yocta::Logger logger;
}

// Project includes
#include "DOT_BasicTypes.h"
#include "DOT_Histogram2D.h"
#include "DOT_Histogram.h"
#include "DOT_Solvers.h"
#include "DOT_CG.h"

using namespace DOT;

// Default test: Wasserstein distance order 1, ground distance L1, two classic images
void single_test(
   const std::string& f1,
   const std::string& f2,
   const GroundDistance& gd,
   bool exact,
   int L
) {
   Histogram2D h1(f1);
   Histogram2D h2(f2);

   logger.info("Total weight H1: %ld, H2: %ld", h1.computeTotalWeight(), h2.computeTotalWeight());

   // Set up configuration option
   Config config;
   config.algo = Algorithm::FlowSimplex;
   config.ground_dist = gd;

   if (gd == GroundDistance::L2) {
      if (exact) {
         int n = static_cast<int>(h1.getN());
         config.buildCoprimes(n - 1);
      } else
         config.buildCoprimes(L);
      logger.info("Size of A_L: %lu", config.coprimes.size());
   }

   // Time vars
   using namespace std;
   using namespace std::chrono;
   // Start time
   time_point<system_clock> start = system_clock::now();

   int64_t wd1 = compute_wd1(h1, h2, config);

   // End time
   time_point<system_clock> end = system_clock::now();
   duration<double> inlineTimeElapsed = end - start;
#ifdef _WIN32
   logger.note("Distance: %I64d, Time: %ld ms, Ground: %s, %s, %s",
               wd1, duration_cast<milliseconds>(inlineTimeElapsed).count(), gd_to_string(config.ground_dist).c_str(),
               f1.c_str(), f2.c_str());
#else
   logger.note("Distance: %ld, Time: %ld ms, Ground: %s, %s, %s",
               wd1, duration_cast<milliseconds>(inlineTimeElapsed).count(), gd_to_string(config.ground_dist).c_str(),
               f1.c_str(), f2.c_str());
#endif
}

// Default test: Wasserstein distance order 1, ground distance L1, two classic images
void single_test_bipartite_W2(
   const std::string& f1,
   const std::string& f2,
   const GroundDistance& gd
) {
   Histogram2D h1(f1);
   Histogram2D h2(f2);

   logger.info("Total weight H1: %ld, H2: %ld", h1.computeTotalWeight(), h2.computeTotalWeight());

   // Set up configuration option
   Config config;
   config.algo = Algorithm::NetworkSimplexTripartite;
   if (gd == GroundDistance::L1)
      config.algo = Algorithm::NetworkSimplexBipartite;

   // Time vars
   using namespace std;
   using namespace std::chrono;
   // Start time
   time_point<system_clock> start = system_clock::now();

   int64_t wd = compute_wd2(h1, h2, config);

   // End time
   time_point<system_clock> end = system_clock::now();
   duration<double> inlineTimeElapsed = end - start;
#ifdef _WIN32
   logger.note("Distance: %I64d, %f, Time: %ld ms, Ground: %s, %s, %s",
               wd, double(wd)/(h1.computeTotalWeight())*32*32/(h1.getN()*h1.getN()), duration_cast<milliseconds>(inlineTimeElapsed).count(),
               gd_to_string(config.ground_dist).c_str(),
               f1.c_str(), f2.c_str());
#else
   logger.note("Distance: %ld, %f, Time: %ld ms, Ground: %s, %s, %s",
               wd, double(wd) / (h1.computeTotalWeight())* 32*32 / (h1.getN()*h1.getN()), duration_cast<milliseconds>(inlineTimeElapsed).count(),
               gd_to_string(config.ground_dist).c_str(),
               f1.c_str(), f2.c_str());
#endif
}

// Default test: Wasserstein distance order 1, ground distance L1, two classic images
void single_test_flowcytometry_W2(
   const std::string& f1,
   const std::string& f2,
   const GroundDistance& gd,
   size_t n,
   size_t d
) {
   Histogram h1(f1, n, d);
   Histogram h2(f2, n, d);

   logger.info("Total weight H1: %ld, H2: %ld", h1.computeTotalWeight(), h2.computeTotalWeight());

   // Set up configuration option
   Config config;
   config.algo = Algorithm::NetworkSimplexTripartite;
   if (gd == GroundDistance::L1)
      config.algo = Algorithm::NetworkSimplexBipartite;

   // Time vars
   using namespace std;
   using namespace std::chrono;
   // Start time
   time_point<system_clock> start = system_clock::now();

   int64_t wd = compute_wd2_fcs(h1, h2, config);

   // End time
   time_point<system_clock> end = system_clock::now();
   duration<double> inlineTimeElapsed = end - start;
#ifdef _WIN32
   logger.note("Distance: %I64d, %f, Time: %ld ms, Ground: %s, %s, %s",
               wd, double(wd) / (h1.computeTotalWeight()) * 32 * 32 / (h1.getN()*h1.getN()), duration_cast<milliseconds>(inlineTimeElapsed).count(),
               gd_to_string(config.ground_dist).c_str(),
               f1.c_str(), f2.c_str());
#else
   logger.note("Distance: %ld, %f, Time: %ld ms, Ground: %s, %s, %s",
               wd, double(wd) / (h1.computeTotalWeight()) * 32 * 32 / (h1.getN()*h1.getN()), duration_cast<milliseconds>(inlineTimeElapsed).count(),
               gd_to_string(config.ground_dist).c_str(),
               f1.c_str(), f2.c_str());
#endif
}

// Default test: Wasserstein distance order 1, ground distance L1, two classic images
void all_test_bipartite_W2() {
   // Base test
   std::string base = "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";
   std::string SEP = "\\";
#ifdef CINECA
   base = "../data/DOTmark/";
   SEP = "/";
#endif

   std::vector<std::string> dirs = { "ClassicImages", "MicroscopyImages", "Shapes" };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = { "32", "64", "128" };
   Config config;

   for (const auto& algo : {
   Algorithm::NetworkSimplexTripartite, Algorithm::NetworkSimplexBipartite
}) {
      config.algo = algo;
      for (const auto& S : Ss) {
         for (const auto& dtype : dirs) {
            for (const auto& f11 : Fs) {
               for (const auto& f22 : Fs)
                  if (f11 < f22) {
                     std::string f1 = "data" + S + "_" + f11;
                     std::string f2 = "data" + S + "_" + f22;

                     Histogram2D h1(base + dtype + SEP + f1);
                     Histogram2D h2(base + dtype + SEP + f2);

                     logger.debug("Total weight H1: %ld, H2: %ld", h1.computeTotalWeight(), h2.computeTotalWeight());

                     // Time vars
                     using namespace std;
                     using namespace std::chrono;
                     // Start time
                     time_point<system_clock> start = system_clock::now();

                     int64_t wd = compute_wd2(h1, h2, config);

                     // End time
                     time_point<system_clock> end = system_clock::now();
                     duration<double> inlineTimeElapsed = end - start;

                     logger.note("Distance: %I64d %f, Time: %ld ms, Algo: %s, %s, %s, %s",
                                 wd, double(wd) / (h1.computeTotalWeight()) * 32 * 32 / (h1.getN()*h1.getN()),
                                 duration_cast<milliseconds>(inlineTimeElapsed).count(), algo_to_string(algo).c_str(),
                                 dtype.c_str(), f1.c_str(), f2.c_str());
                  }
            }
         }
      }
   }
}


// Default test: Wasserstein distance order 1, ground distance L1, two classic images
void default_test() {
   // Base test
   std::string base = "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";
   std::string SEP = "\\";
#ifdef CINECA
   base = "../data/DOTmark/";
   SEP = "/";
#endif

   std::string dtype = "Shapes";

   std::string f1 = "data64_1004.csv";
   std::string f2 = "data64_1008.csv";

   Histogram2D h1(base + dtype + SEP + f1);
   Histogram2D h2(base + dtype + SEP + f2);

   logger.info("Total weight H1: %ld, H2: %ld", h1.computeTotalWeight(), h2.computeTotalWeight());

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
   using namespace std;
   using namespace std::chrono;
   // Start time
   time_point<system_clock> start = system_clock::now();

   int64_t wd1 = compute_wd1(h1, h2, config);

   // End time
   time_point<system_clock> end = system_clock::now();
   duration<double> inlineTimeElapsed = end - start;

#ifdef _WIN32
   logger.note("Distance: %I64d, Time: %ld ms, Ground: %s, %s, %s",
               wd1, duration_cast<milliseconds>(inlineTimeElapsed).count(), gd_to_string(config.ground_dist).c_str(),
               f1.c_str(), f2.c_str());
#else
   logger.note("Distance: %ld, Time: %ld ms, Ground: %s, %s, %s",
               wd1, duration_cast<milliseconds>(inlineTimeElapsed).count(), gd_to_string(config.ground_dist).c_str(),
               f1.c_str(), f2.c_str());
#endif

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

   std::vector<std::string> Ss = { "32" };//, "64", "128", "256", "512" };

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
//         } else {
//            config.buildCoprimes(2);
         }
         for (const auto& dtype : dirs) {
            for (const auto& f11 : Fs) {
               for (const auto& f22 : Fs)
                  if (f11 < f22) {
                     std::string f1 = "data" + S + "_" + f11;
                     std::string f2 = "data" + S + "_" + f22;

                     Histogram2D h1(base + dtype + SEP + f1);
                     Histogram2D h2(base + dtype + SEP + f2);

                     logger.debug("Total weight H1: %ld, H2: %ld", h1.computeTotalWeight(), h2.computeTotalWeight());

                     // Time vars
                     using namespace std;
                     using namespace std::chrono;
                     // Start time
                     time_point<system_clock> start = system_clock::now();

                     int64_t wd1 = compute_wd1(h1, h2, config);

                     // End time
                     time_point<system_clock> end = system_clock::now();
                     duration<double> inlineTimeElapsed = end - start;

                     logger.note("Distance: %I64d, Time: %ld ms, Ground: %s, %s, %s, %s",
                                 wd1, duration_cast<milliseconds>(inlineTimeElapsed).count(), gd_to_string(config.ground_dist).c_str(),
                                 dtype.c_str(), f1.c_str(), f2.c_str());
                  }
            }
         }
      }
   }
}

// Default test: Wasserstein distance order 1, ground distance L1, two classic images
void all_fcs_test() {
   // Base test
   std::string base = "D:\\Ricerca\\DOTA\\data\\AML\\CSV\\";
   std::string SEP = "\\";

   std::vector<std::string> Fs = { "0001.csv", "0009.csv", "0017.csv", "0025.csv", "0033.csv",
                                   "0041.csv", "0049.csv", "0057.csv", "0065.csv", "0201.csv"
                                 };

   std::vector<Algorithm> Gs = { DOT::Algorithm::NetworkSimplexBipartite };

   std::vector<size_t> Ds = { 2, 3 };

   int n = 32;

   // Set up configuration option
   Config config;

   for (Algorithm gd : Gs) {
      config.algo = gd;
      for (const auto& d : Ds)
         if (!(d == 4 && gd == Algorithm::NetworkSimplexBipartite)) {
            for (const auto& f1 : Fs) {
               for (const auto& f2 : Fs)
                  if (f1 < f2) {
                     Histogram h1(base + SEP + f1, n, d);
                     Histogram h2(base + SEP + f2, n, d);

                     logger.debug("Total weight H1: %ld, H2: %ld", h1.computeTotalWeight(), h2.computeTotalWeight());

                     // Time vars
                     using namespace std;
                     using namespace std::chrono;
                     // Start time
                     time_point<system_clock> start = system_clock::now();

                     int64_t wd = compute_wd2_fcs(h1, h2, config);

                     // End time
                     time_point<system_clock> end = system_clock::now();
                     duration<double> inlineTimeElapsed = end - start;

                     logger.note("Distance: %I64d, Time: %ld ms, n: %d, d: %d, Algo: %s, %s, %s",
                                 wd, duration_cast<milliseconds>(inlineTimeElapsed).count(), int(n), int(d),
                                 algo_to_string(config.algo).c_str(),
                                 f1.c_str(), f2.c_str());
                  }
            }
         }
   }
}

// Column Generation test
void column_generation_test() {
   int n = 128;
   int seed = 13;
   auto X = uniformSpace2D(n, ++seed);
   auto Y = uniformSpace2D(n, ++seed);

   auto PI = randomMarginal(n*n, 255, ++seed);
   auto MU = randomMarginal(n*n, 255, ++seed);

   ColumnGeneration(X, PI, Y, MU);
}


// Main entry point
int main(int argc, char* argv[]) {
   args::ArgumentParser parser("DOTLib v0.3.0: Discrete Optimal Transport library.", "");
   args::HelpFlag help(parser, "help", "Display this help menu", { 'h', "help" });
   args::ValueFlag<std::string> h1_filename(parser, "h1_filename", "filename of the first histogram", { 'f', "h1" });
   args::ValueFlag<std::string> h2_filename(parser, "h2_filename", "filename of the second histogram", { 'g', "h2" });
   args::ValueFlag<std::string> log_filename(parser, "log_filename", "logger filename", { 'l', "log" });

   args::ValueFlag<int> sides(parser, "sides", "Dimension of the subgrid size L", { 'L', "sides" });
   args::ValueFlag<int> verbosity(parser, "verbosity", "Verbosity level of logger", { 'v', "ver" });

   args::Flag dotmarks(parser, "dotmarks", "Run massive tests on DotMarks instances", { 'd', "dotmarks" });
   args::Flag single(parser, "single", "Run single test between two histograms", { 's', "single" });
   args::Flag exact(parser, "exact", "When using the L2 norm as ground distance, solve the exact problem", { 'e', "exact" });
   args::Flag colgen(parser, "columnGeneration", "Solve the problem by column generation", { 'z', "colgen" });

   args::Group nomrs(parser, "This group is all exclusive:", args::Group::Validators::AtMostOne);
   args::Flag norm1(nomrs, "L1", "Norm-1 ground distance", { 'a', "l1" });
   args::Flag norm2(nomrs, "L2", "Norm-2 ground distance", { 'b', "l2" });
   args::Flag norm8(nomrs, "Linf", "Norm-Inf ground distance", { 'c', "l8" });

   try {
      if (argc == 1) {
         std::cout << argv[0] << parser;
         return EXIT_SUCCESS;
      }

      parser.ParseCLI(argc, argv);

   } catch (args::Help) {
      std::cout << parser;
      return EXIT_SUCCESS;
   } catch (args::ParseError e) {
      std::cerr << e.what() << std::endl;
      std::cerr << parser;
      return EXIT_FAILURE;
   } catch (args::ValidationError e) {
      std::cerr << e.what() << std::endl;
      std::cerr << parser;
      return EXIT_FAILURE;
   }

   // Start the command line options
   try {
      // Set logger first
      if (log_filename)
         logger.setFileStream(args::get(log_filename));
      if (verbosity)
         logger.setVerbosityLevel(static_cast<yocta::VerbosityLevel>(args::get(verbosity)));

      // Run massive test on DOTMARKS data set
      if (dotmarks) {
         all_dotmark_test();
         exit(EXIT_SUCCESS);
      }

      // Run Column Generation Algorithm
      if (colgen) {
         column_generation_test();
         exit(EXIT_SUCCESS);
      }

      // Run single test
      if (single) {
         std::string f1 = "data32_1001.csv";
         std::string f2 = "data32_1002.csv";
         if (h1_filename && h2_filename) {
            f1 = args::get(h1_filename);
            f2 = args::get(h2_filename);
         }
         GroundDistance gd = GroundDistance::L1;
         if (norm2)
            gd = GroundDistance::L2;
         if (norm8)
            gd = GroundDistance::Linf;

         // Run single test
         single_test(f1, f2, gd, exact, 1);

         // FCS - NIPS 2018
         //single_test_flowcytometry_W2(f1,f2,gd);
         //all_fcs_test();

         exit(EXIT_SUCCESS);
      }
   } catch (std::exception e) {
      logger.error("Runtime error: %s", e.what());
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
