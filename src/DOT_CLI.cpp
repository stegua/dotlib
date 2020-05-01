/*
*  Main authors:
*     @author Stefano Gualandi <stefano.gualandi@gmail.com>
*
*     @fileoverview Copyright (c) 2017-19, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
*  Last update: May, 2019
*/


// Yocta Logger library
#include <yocta_logger.hh>

// ----------------- Define ALL TEST CASES -------------------
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_DISABLE_STRINGIFICATION
#include "catch.hpp"

#include "DOT_BipartiteGurobi.hpp"
#include "DOT_BipartiteLemon.hpp"
#include "DOT_KWD_P1.hpp"


// Base test
std::string SEP = "\\";

std::string base = "C:\\Users\\Gualandi\\Google Drive\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";
//std::string base = "D:\\GoogleDrive\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";


TEST_CASE("BIPARTITE_GUROBI") {
   std::vector<std::string> dirs = { "ClassicImages" };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", //"1003.csv", "1004.csv", "1005.csv",
                                   //"1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32", "64"
   };

   for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                  DOT::BipartiteLP(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " GRB_PRIMAL");
                  DOT::BipartiteLP(Mu, Nu, 1, dtype + " " + f11 + " " + f22 + " " + S + " GRB_DUAL");
                  DOT::BipartiteLP(Mu, Nu, 2, dtype + " " + f11 + " " + f22 + " " + S + " GRB_BARRIER");
               }
         }
      }
   }
}


TEST_CASE("BIPARTITE_LEMON_SIMPLEX") {
   std::vector<std::string> dirs = { "ClassicImages", //"WhiteNoise", "CauchyDensity", "GRFmoderate", "GRFrough", "GRFsmooth",
                                     //"LogGRF", "LogitGRF", "MicroscopyImages", "Shapes",
                                   };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "64"
   };

   for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                  DOT::BipartiteNetworkSimplex(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " LEMON_NS");
                  DOT::BipartiteNetworkSimplex(Mu, Nu, 1, dtype + " " + f11 + " " + f22 + " " + S + " LEMON_NS");
                  DOT::BipartiteNetworkSimplex(Mu, Nu, 2, dtype + " " + f11 + " " + f22 + " " + S + " LEMON_NS");
               }
         }
      }
   }
}

TEST_CASE("BIPARTITE_LEMON_COST_SCALING") {
   std::vector<std::string> dirs = { "ClassicImages", //"WhiteNoise", "CauchyDensity", "GRFmoderate", "GRFrough", "GRFsmooth",
                                     //"LogGRF", "LogitGRF", "MicroscopyImages", "Shapes",
                                   };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", //"1003.csv", //"1004.csv", "1005.csv",
                                   //"1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32", "64"
   };

   for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                  DOT::BipartiteCostScaling(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " LEMON_CS");
               }
         }
      }
   }
}

TEST_CASE("BIPARTITE_LEMON_CYCLE_CANCELING") {
   std::vector<std::string> dirs = { "ClassicImages" };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32", "64"
   };

   for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                  DOT::BipartiteCycleCanceling(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " LEMON_CC");
               }
         }
      }
   }
}

TEST_CASE("BIPARTITE_EMD") {
   std::vector<std::string> dirs = { "ClassicImages" };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32"
   };

   for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                  DOT::BipartiteEMD(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " EMD");
               }
         }
      }
   }
}

TEST_CASE("EMD_L1_LINK_OKADA") {
   std::vector<std::string> dirs = { "ClassicImages" };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", //"1003.csv", "1004.csv", "1005.csv",
                                   //"1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32", "64", "128", "256", "512"
   };

   for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                  DOT::EMD_L1(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " EMD_L1");
               }
         }
      }
   }
}


TEST_CASE("BGV_L1") {
   std::vector<std::string> dirs = { "ClassicImages", //"WhiteNoise", "CauchyDensity", "GRFmoderate", "GRFrough", "GRFsmooth",
                                     //"LogGRF", "LogitGRF", "MicroscopyImages", "Shapes",
                                   };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", //"1003.csv", "1004.csv", "1005.csv",
                                   //"1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32", "64", "128", "256", "512"
   };

   for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                  DOT::BGV_L1(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " BGV_L1");
                  DOT::BGV_L1(Mu, Nu, 1, dtype + " " + f11 + " " + f22 + " " + S + " BGV_L1_cs");
               }
         }
      }
   }
}

TEST_CASE("BGV_Linf") {
   std::vector<std::string> dirs = { "ClassicImages", //"WhiteNoise", "CauchyDensity", "GRFmoderate", "GRFrough", "GRFsmooth",
                                     //"LogGRF", "LogitGRF", "MicroscopyImages", "Shapes",
                                   };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", //"1003.csv", "1004.csv", "1005.csv",
                                   //"1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32", "64", "128", "256", "512"
   };

   for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                  DOT::BGV_Linf(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " BGV_Linf");
                  DOT::BGV_Linf(Mu, Nu, 1, dtype + " " + f11 + " " + f22 + " " + S + " BGV_Linf_cs ");
               }
         }
      }
   }
}

TEST_CASE("BGV_L2") {
   std::vector<std::string> dirs = { "ClassicImages", "WhiteNoise", "CauchyDensity", "GRFmoderate", "GRFrough", "GRFsmooth",
                                     "LogGRF", "LogitGRF", "MicroscopyImages", "Shapes",
                                   };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      //"32", "64",
      "128",
   };

   for (const auto& S : Ss) {
      size_t s = static_cast<size_t>(stoi(S));

      auto coprimes = DOT::buildCoprimes(s-1);

      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                  DOT::BGV_L2(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " BGV_L2", coprimes);
               }
         }
      }
   }
}

TEST_CASE("BGV_L2_APPROX") {
   std::vector<std::string> dirs = { "ClassicImages", "WhiteNoise", "CauchyDensity", "GRFmoderate", "GRFrough", "GRFsmooth",
                                     "LogGRF", "LogitGRF", "MicroscopyImages", "Shapes",
                                   };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32", "64", "128", "256", "512"
   };

   for (int L : {
            2, 3, 5, 10
         }) {
      // Build coprime subsets
      auto coprimes = DOT::buildCoprimes(L);

      for (const auto& S : Ss) {
         for (const auto& dtype : dirs) {
            for (const auto& f11 : Fs) {
               for (const auto& f22 : Fs)
                  if (f11 < f22) {
                     std::string f1 = "data" + S + "_" + f11;
                     std::string f2 = "data" + S + "_" + f22;

                     DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                     DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                     DOT::BGV_L2(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " BGV_L2 " + std::to_string(L), coprimes);
                  }
            }
         }
      }
   }
}


//-------------------------------------------------------------------
// TEST FOR PREPARING THE TABLE FOR THE PAPER
//-------------------------------------------------------------------

TEST_CASE("SIAMOPT_TAB1") {
   std::vector<std::string> dirs = { "ClassicImages" };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32"
   };

   for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                  DOT::BipartiteEMD(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " EMD");

                  DOT::BipartiteLP(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " GRB_PRIMAL");
                  DOT::BipartiteLP(Mu, Nu, 1, dtype + " " + f11 + " " + f22 + " " + S + " GRB_DUAL");
                  DOT::BipartiteLP(Mu, Nu, 2, dtype + " " + f11 + " " + f22 + " " + S + " GRB_BARRIER");

                  DOT::BipartiteCycleCanceling(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " LEMON_CC");
                  DOT::BipartiteCostScaling(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " LEMON_CS");
                  DOT::BipartiteNetworkSimplex(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " LEMON_NS");
               }
         }
      }
   }
}

TEST_CASE("SIAMOPT_TAB2") {
   std::vector<std::string> dirs = { "ClassicImages" };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32","64", "128", "256", "512"
   };

   for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);


                  DOT::BGV_L1(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " BGV_L1");
                  DOT::BGV_L1(Mu, Nu, 1, dtype + " " + f11 + " " + f22 + " " + S + " BGV_L1_cs");
                  DOT::EMD_L1(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " EMD_L1");
               }
         }
      }
   }
}

TEST_CASE("SIAMOPT_TAB3") {
   std::vector<std::string> dirs = { "ClassicImages", "WhiteNoise", "CauchyDensity", "GRFmoderate", "GRFrough", "GRFsmooth",
                                     "LogGRF", "LogitGRF", "MicroscopyImages", "Shapes"
                                   };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32","64", "128", "256", "512"
   };

   for (const auto& S : Ss) {
      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);


                  DOT::BGV_L1(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " BGV_L1");
                  DOT::BGV_Linf(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " BGV_Linf");
               }
         }
      }
   }

   std::vector<std::string> S2= {
      "32", "64", "128",
   };

   for (const auto& S : S2) {
      size_t s = static_cast<size_t>(stoi(S));

      auto coprimes = DOT::buildCoprimes(s - 1);

      for (const auto& dtype : dirs) {
         for (const auto& f11 : Fs) {
            for (const auto& f22 : Fs)
               if (f11 < f22) {
                  std::string f1 = "data" + S + "_" + f11;
                  std::string f2 = "data" + S + "_" + f22;

                  DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                  DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                  DOT::BGV_L2(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " BGV_L2", coprimes);
               }
         }
      }
   }
}

TEST_CASE("SIAMOPT_TAB4") {
   std::vector<std::string> dirs = { "ClassicImages", "WhiteNoise", "CauchyDensity", "GRFmoderate", "GRFrough", "GRFsmooth",
                                     "LogGRF", "LogitGRF", "MicroscopyImages", "Shapes",
                                   };

   std::vector<std::string> Fs = { "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
                                   "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"
                                 };

   std::vector<std::string> Ss = {
      "32", "64", "128", "256", "512"
   };

   for (const auto& S : Ss) {
      for (int L : {
               2, 3, 5, 10
            }) {
         // Build coprime subsets
         auto coprimes = DOT::buildCoprimes(L);

         for (const auto& dtype : dirs) {
            for (const auto& f11 : Fs) {
               for (const auto& f22 : Fs)
                  if (f11 < f22) {
                     std::string f1 = "data" + S + "_" + f11;
                     std::string f2 = "data" + S + "_" + f22;

                     DOT::MeasureR2 Mu(base + dtype + SEP + f1);
                     DOT::MeasureR2 Nu(base + dtype + SEP + f2);

                     DOT::BGV_L2(Mu, Nu, 0, dtype + " " + f11 + " " + f22 + " " + S + " BGV_L2 " + std::to_string(L), coprimes);
                  }
            }
         }
      }
   }
}