#include "OT_Solvers.h""

#include "OT_Exact_L2.h"
#include "OT_Utils.h"

// Test for a 3x3 binary matrix
real_t test_3_3(int norm) {
   histogram_t image1 = { 1, 0, 0, 0, 0, 0, 0, 0, 0 };
   histogram_t image2 = { 0, 0, 0, 0, 0, 0, 0, 1, 0 };

   real_t cost = std::numeric_limits<real_t>::max();

   if (norm == 1) {
      cost = DOTLib::solve_L1_1(image1, image2);
      spd::get("console")->info("Optimal cost L_1_1: {:03.3f}", cost);
   }
   if (norm == 1000) {
      cost = DOTLib::solve_Linf_1(image1, image2);
      spd::get("console")->info("Optimal cost L_inf_1: {:03.3f}", cost);
   }

   return cost;
}

// Test for
real_t test_W_1(const std::string& n_str, int d) {
   size_t n = std::stoi(n_str.c_str());

   histogram_t image1 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data" + n_str + "_1005.csv");
   histogram_t image2 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data" + n_str + "_1009.csv");

   // Time vars
   std::chrono::time_point<std::chrono::system_clock> start, end;
   // Start time.
   start = std::chrono::system_clock::now();
   real_t exact_cost = 0;

   if (d == 1)
      exact_cost = solve_exact_W_1_1(image1, image2);
   if (d == 2)
      exact_cost = solve_exact_W_1_2(image1, image2);
   if (d == 8)
      exact_cost = solve_exact_W_1_8(image1, image2);

   // End time.
   end = std::chrono::system_clock::now();
   std::chrono::duration<double> inlineTimeElapsed = end - start;

   // End time.
   end = std::chrono::system_clock::now();
   inlineTimeElapsed = end - start;
   spd::get("console")->info("Cost: {:.17f} - Time: {} ms", exact_cost,  std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count());

   return exact_cost;
}

// Test for
real_t test_all_exact_d2(const std::string& n_str) {
   size_t n = std::stoi(n_str.c_str());

   std::string base = "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data";
   //std::string base = "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\LogGRF\\";

   std::vector<std::string> Fs = { "_1001.csv", "_1002.csv", "_1003.csv", "_1004.csv", "_1005.csv",
                                   "_1006.csv", "_1007.csv", "_1008.csv", "_1009.csv", "_1010.csv"
                                 };

   for (const auto& e1 : Fs)
      for (const auto& e2 : Fs) {
         if (e1 < e2) {
            auto f1 = n_str + e1;
            auto f2 = n_str + e2;

            histogram_t image1 = read_image(n, base + f1);
            histogram_t image2 = read_image(n, base + f2);

            spd::get("console")->info("{} {}", f1, f2);

            // Time vars
            std::chrono::time_point<std::chrono::system_clock> start, end;
            // Start time.
            start = std::chrono::system_clock::now();
            real_t exact_cost = 0;

            exact_cost = solve_exact_W_1_2(image1, image2);

            // End time.
            end = std::chrono::system_clock::now();
            std::chrono::duration<double> inlineTimeElapsed = end - start;

            // End time.
            end = std::chrono::system_clock::now();
            inlineTimeElapsed = end - start;
            spd::get("console")->info("Cost: {:.17f} - Time: {} ms", exact_cost, std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count());

         }
      }

   return 1;
}

real_t test_all_heuristic_d2(const std::string& n_str) {
   size_t n = std::stoi(n_str.c_str());

   std::string base = "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";

   std::vector<std::string> dirs = { "CauchyDensity", "ClassicImages", "GRFmoderate", "GRFrough", "GRFsmooth",
                                     "LogGRF", "LogitGRF", "MicroscopyImages", "Shapes", "WhiteNoise"
                                   };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<int> Ls = {2,3,5,10};
   Ls.push_back(n - 1);

   for (const std::string& ds : dirs) {
      for (const auto& e1 : Fs)
         for (const auto& e2 : Fs)
            if (e1 < e2) {
               for (int L : Ls) {
                  auto f1 = n_str + "_" + e1;
                  auto f2 = n_str + "_" + e2;

                  spd::get("logs")->info("{} {} {}", ds, f1, f2);

                  histogram_t image1 = read_image(n, base + ds + "\\data" + f1);
                  histogram_t image2 = read_image(n, base + ds + "\\data" + f2);

                  // Time vars
                  std::chrono::time_point<std::chrono::system_clock> start, end;
                  // Start time.
                  start = std::chrono::system_clock::now();
                  real_t cost = 0;

                  cost = DOTLib::solve_L2_1(image1, image2, L);

                  // End time.
                  end = std::chrono::system_clock::now();
                  std::chrono::duration<double> inlineTimeElapsed = end - start;

                  // End time.
                  end = std::chrono::system_clock::now();
                  inlineTimeElapsed = end - start;
                  spd::get("logs")->info("Cost: {:.17f} - Time: {} ms", cost, std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count());

               }
               break;
            }
   }

   return 1;
}

// Test for
real_t test_V_1(const std::string& n_str, int d) {
   size_t n = std::stoi(n_str.c_str());

   auto f1 = n_str + "_1005";
   auto f2 = n_str + "_1009";

   histogram_t image1 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data" + f1 + ".csv");
   histogram_t image2 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data" + f2 + ".csv");

   spd::get("console")->info("{} {}", f1, f2);

   // Time vars
   std::chrono::time_point<std::chrono::system_clock> start, end;
   // Start time.
   start = std::chrono::system_clock::now();
   real_t exact_cost = 0;

   if (d == 1)
      exact_cost = DOTLib::solve_L1_1(image1, image2);
   if (d == 2)
      exact_cost = DOTLib::solve_L2_1(image1, image2, n-1);
   if (d == 8)
      exact_cost = DOTLib::solve_Linf_1(image1, image2);


   // End time.
   end = std::chrono::system_clock::now();
   std::chrono::duration<double> inlineTimeElapsed = end - start;

   // End time.
   end = std::chrono::system_clock::now();
   inlineTimeElapsed = end - start;
   spd::get("console")->info("Cost: {:.2f} - Time: {} ms", exact_cost, std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count());

   return exact_cost;
}

// ----------------- Define ALL TEST CASES -------------------
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

TEST_CASE("Optimal transport are solved", "[EMD_p_1]") {
   // Console logger with color
   auto console = spd::stdout_color_mt("console");
   console->info("DOTlib v0.1.0");
   //if (false) {
   //   // Test L_1_1 norm on a 3x3 matrix
   //   REQUIRE(fabs(test_3_3(1) - 3.0) < 1e-09);

   //   // Test L_infinity_1 norm on a 3x3 matrix
   //   REQUIRE(fabs(test_3_3(1000) - 2.0) < 1e-09);

   //   // Test on bipartite graph

   //   REQUIRE(test_W_1("32", 2) >= 238308025);
   //   REQUIRE(test_W_1("64", 2) >= 1906840961);

   //   REQUIRE(test_W_1("32", 1) >= 238308025);
   //   REQUIRE(test_W_1("64", 1) >= 1906840961);

   //   REQUIRE(test_W_1("32", 8) >= 0);
   //   REQUIRE(test_W_1("64", 8) >= 0);

   //   REQUIRE(test_V_1("32", 1) >= 0);
   //   REQUIRE(test_V_1("64", 1) >= 0);
   //   REQUIRE(test_V_1("128", 1) >= 0);
   //   //REQUIRE(test_V_1("256", 1) >= 0);
   //   //REQUIRE(test_V_1("512", 1) >= 0);

   //   REQUIRE(test_V_1("32", 8) >= 0);
   //   REQUIRE(test_V_1("64", 8) >= 0);
   //   REQUIRE(test_V_1("128", 8) >= 0);
   //   //REQUIRE(test_V_1("256", 8) >= 0);
   //   //REQUIRE(test_V_1("512", 8) >= 0);
   //}


   //REQUIRE(test_W_1("64", 2) >= 238308025);
   /*REQUIRE(test_V_1("32", 2) >= 0);
   REQUIRE(test_V_1("64", 2) >= 0);
   REQUIRE(test_V_1("128", 2) >= 0);
   REQUIRE(test_V_1("256", 2) >= 0);
   REQUIRE(test_V_1("512", 2) >= 0);*/

   // All pair exact
   //REQUIRE(test_all_exact_d2("32") >= 0);
   // REQUIRE(test_all_exact_d2("64") >= 0);

   auto my_logger = spd::basic_logger_mt("logs", "./log_test_32_64_128_all_h.csv");

   REQUIRE(test_all_heuristic_d2("32") >= 0);
   REQUIRE(test_all_heuristic_d2("64") >= 0);
   REQUIRE(test_all_heuristic_d2("128") >= 0);

}
