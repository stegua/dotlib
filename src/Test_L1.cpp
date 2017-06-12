#include "OT_Solvers.h""

// Test for a 3x3 binary matrix
real_t test_3_3() {
   histogram_t image1 = { 1, 0, 0, 0, 0, 0, 0, 0, 0 };
   histogram_t image2 = { 0, 0, 0, 0, 0, 0, 0, 1, 0 };

   real_t cost = solve_L1_1(image1, image2);

   spd::get("console")->info("Optimal cost L1_1: {:03.3f}", cost);

   return cost;
}


// ----------------- Define ALL TEST CASES -------------------
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

TEST_CASE("Factorials are computed", "[factorial]") {
   // Console logger with color
   auto console = spd::stdout_color_mt("console");
   console->info("DOTlib v0.0.1");

   REQUIRE(fabs(test_3_3() - 2.0) < 1e-09);
}