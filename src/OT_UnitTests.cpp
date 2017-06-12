#include "OT_Solvers.h""

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


// ----------------- Define ALL TEST CASES -------------------
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

TEST_CASE("Optimal transport are solved", "[EMD_p_1]") {
   // Console logger with color
   auto console = spd::stdout_color_mt("console");
   console->info("DOTlib v0.0.1");

   // Test L_1_1 norm on a 3x3 matrix
   REQUIRE(fabs(test_3_3(1) - 3.0) < 1e-09);

   // Test L_infinity_1 norm on a 3x3 matrix
   REQUIRE(fabs(test_3_3(1000) - 2.0) < 1e-09);
}