/**
* @fileoverview Copyright (c) 2017, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#include <chrono>

#include "OT_BasicTypes.h"
#include "OT_Apx_L2_1.h"
#include "OT_Exact_L2.h"
#include "OT_Utils.h"

/**
* @brief Read images from text file and return corresponfing matrix
*/
matrix_t compute_distance(size_t n, int p = 2, int q = 1) {
   matrix_t cost;

   auto m = [&n](size_t x, size_t y) {
      return x * n + y;
   };

   for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
         cost.emplace_back(row_t(n * n, 0));
         for (int v = 0; v < n; ++v)
            for (int w = 0; w < n; ++w)
               cost[m(i, j)][m(v, w)] = std::trunc(65536 * pow(pow(pow(fabs(i - v), p) + pow(fabs(j - w), p), 1.0 / p), q)) / 65536;
      }

   return cost;
}

/**
* @brief Compute Earth Moving Distance (EMD) between pair of images
*/
int compute_EMD(const histogram_t& h1, const histogram_t& h2, const matrix_t& cost) {
   using namespace lemon;

   real_t distance = std::numeric_limits<real_t>::max();

   size_t d = h1.size();

   // Build the graph for max flow
   Graph g;

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<Graph::Node> nodes;
   // add first d source nodes
   for (size_t i = 0; i < d; ++i)
      nodes.emplace_back(g.addNode());
   // add first d destination nodes
   for (size_t i = 0; i < d; ++i)
      nodes.emplace_back(g.addNode());

   std::vector<Graph::Arc> arcs;
   for (size_t i = 0; i < d; ++i)
      for (size_t j = 0; j < d; ++j) {
         Graph::Arc a = g.addArc(nodes[i], nodes[d + j]);
         arcs.emplace_back(a);
      }

   // graph must be finished before initializing cycle canceling !
   CycleCanceling<Graph, LimitValueType, LimitValueType> cycle(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g), c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i)
      b_i[nodes[i]] = int(h1[i]);
   for (size_t i = 0; i < d; ++i)
      b_i[nodes[d + i]] = -int(h2[i]);

   // Add all edges
   for (size_t i = 0; i < d; ++i)
      for (size_t j = 0; j < d; ++j) {
         const auto& a = arcs[i * d + j];
         l_i[a] = 0;
         u_i[a] = std::min(h1[i], h2[i]);// cycle.INF;
         c_i[a] = static_cast<LimitValueType>(cost[i][j]);  // cost
      }

   //set lower/upper bounds, cost
   cycle.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   fprintf(stdout, "Input graph created with %d nodes and %d arcs\n\n",
           countNodes(g), countArcs(g));

   CycleCanceling<Graph, LimitValueType, LimitValueType>::ProblemType ret = cycle.run();

   switch (ret) {
   case CycleCanceling<Graph>::INFEASIBLE:
      std::cerr << "INFEASIBLE" << std::endl;
      break;
   case CycleCanceling<Graph>::OPTIMAL:
      std::cerr << "OPTIMAL" << std::endl;
      break;
   case CycleCanceling<Graph>::UNBOUNDED:
      std::cerr << "UNBOUNDED" << std::endl;
      break;
   }

   int sol_value = cycle.totalCost();

   return sol_value;
}

/**
* @brief Compute Earth Moving Distance (EMD) between pair of images
*/
real_t compute_ns_EMD(const histogram_t& h1, const histogram_t& h2, const matrix_t& cost) {
   using namespace lemon;

   real_t distance = std::numeric_limits<real_t>::max();

   size_t d = h1.size();

   // Build the graph for max flow
   Graph g;

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<Graph::Node> nodes;
   // add first d source nodes
   for (size_t i = 0; i < d; ++i)
      nodes.emplace_back(g.addNode());
   // add first d destination nodes
   for (size_t i = 0; i < d; ++i)
      nodes.emplace_back(g.addNode());

   std::vector<Graph::Arc> arcs;
   arcs.reserve(d * d);
   for (size_t i = 0; i < d; ++i)
      for (size_t j = 0; j < d; ++j) {
         Graph::Arc a = g.addArc(nodes[i], nodes[d + j]);
         arcs.emplace_back(a);
      }

   NetworkSimplex<Graph, LimitValueType, real_t> cycle(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<real_t> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i)
      b_i[nodes[i]] = int(h1[i]);
   for (size_t i = 0; i < d; ++i)
      b_i[nodes[d + i]] = -int(h2[i]);

   // Add all edges
   for (size_t i = 0; i < d; ++i)
      for (size_t j = 0; j < d; ++j) {
         const auto& a = arcs[i * d + j];
         l_i[a] = 0;
         u_i[a] = std::min(h1[i], h2[i]);// cycle.INF;
         c_i[a] = static_cast<real_t>(cost[i][j]);  // cost
      }

   //set lower/upper bounds, cost
   cycle.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   fprintf(stdout, "Input graph created with %d nodes and %d arcs\n\n",
           countNodes(g), countArcs(g));

   NetworkSimplex<Graph, LimitValueType, real_t>::ProblemType ret = cycle.run();

   switch (ret) {
   case NetworkSimplex<Graph>::INFEASIBLE:
      std::cerr << "INFEASIBLE" << std::endl;
      break;
   case NetworkSimplex<Graph>::OPTIMAL:
      std::cerr << "OPTIMAL" << std::endl;
      break;
   case NetworkSimplex<Graph>::UNBOUNDED:
      std::cerr << "UNBOUNDED" << std::endl;
      break;
   }

   real_t sol_value = cycle.totalCost();

   return sol_value;
}

/**
* @brief Maint Entry point
*/
int mmmain(int argc, char* argv) {
   size_t n = 32;

   histogram_t image1 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data32_1005.csv");
   histogram_t image2 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data32_1009.csv");

   // Time vars
   std::chrono::time_point<std::chrono::system_clock> start, end;
   // Start time.
   start = std::chrono::system_clock::now();

   auto exact_cost = solve_exact_W_1_2(image1, image2);

   // End time.
   end = std::chrono::system_clock::now();
   std::chrono::duration<double> inlineTimeElapsed = end - start;


   // End time.
   end = std::chrono::system_clock::now();
   inlineTimeElapsed = end - start;
   std::cout << "Cost: " << exact_cost << " - Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count() << " ms \n";

   //std::cout << "Ratio: " << fabs(l21_cost - emd_cost) / emd_cost * 100 << " %\n";

   return 0;
}

#ifdef MY_TEST
int main(int argc, char* argv) {
   // Console logger with color
   //auto console = spd::stdout_color_mt("console");
   //console->info("DOTlib v0.0.1");

   // Formatting examples
   //console->error("Some error message with arg{}..", 1);
   //console->warn("Easy padding in numbers like {:08d}", 12);
   //console->critical("Support for int: {0:d};  hex: {0:x};  oct: {0:o}; bin: {0:b}", 42);
   //console->info("Support for floats {:03.2f}", 1.23456);
   //console->info("Positional args are {1} {0}..", "too", "supported");
   //console->info("{:<30}", "left aligned");


   if (argc == 2) {
      std::string base = "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\";
      //std::string base = "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\LogGRF\\";
      size_t n = 32;
      std::vector<std::string> Fs = { "data32_1001.csv", "data32_1002.csv", "data32_1003.csv", "data32_1004.csv", "data32_1005.csv",
                                      "data32_1006.csv", "data32_1007.csv", "data32_1008.csv", "data32_1009.csv", "data32_1010.csv"
                                    };
      //size_t n = 64;
      //std::vector<std::string> Fs = { "data64_1001.csv", "data64_1002.csv", "data64_1003.csv", "data64_1004.csv", "data64_1005.csv",
      //                                "data64_1006.csv", "data64_1007.csv", "data64_1008.csv", "data64_1009.csv", "data64_1010.csv"
      //                              };

      for (const auto& f1 : Fs)
         for (const auto& f2 : Fs) {
            if (f1 < f2) {
               std::cout << "Comparing: " << (base + f1) << " and " << (base + f2) << std::endl;

               histogram_t image1 = read_image(n, base + f1);
               histogram_t image2 = read_image(n, base + f2);
               // Time vars
               std::chrono::time_point<std::chrono::system_clock> start, end;
               // Start time.
               start = std::chrono::system_clock::now();

               matrix_t cost = compute_distance(n, 2, 1);
               real_t emd_cost = compute_ns_EMD(image1, image2, cost);

               end = std::chrono::system_clock::now();
               std::chrono::duration<double> inlineTimeElapsed = end - start;
               std::cout << "EMD Cost: " << emd_cost << " - Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count() << " ms \n";

               start = std::chrono::system_clock::now();


               auto l21_cost = apx_L2_1(image1, image2);

               // End time.
               end = std::chrono::system_clock::now();
               std::chrono::duration<double> t2 = end - start;
               std::cout << "APX Cost: " << l21_cost << " - Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2).count() << " ms \n";

               std::cout << f1 << " " << f2 << " Ratio: " << fabs(l21_cost - emd_cost) / emd_cost * 100 << " %\n";
            }
         }

      exit(0);
   }

// Size of the square image
   size_t n = 32;
//   histogram_t image1 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data32_1001.csv");
//histogram_t image2 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data32_1005.csv");

   histogram_t image1 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data32_1004.csv");
   histogram_t image2 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data32_1008.csv");

   if (false) {
      n = 3;
      image1 = { 1, 0, 0, 0, 0, 0, 0, 0, 0 };
      image2 = { 0, 0, 0, 0, 0, 0, 0, 1, 0 };
   }

// Time vars
   std::chrono::time_point<std::chrono::system_clock> start, end;
// Start time.
   start = std::chrono::system_clock::now();


// End time.
   end = std::chrono::system_clock::now();
   std::chrono::duration<double> inlineTimeElapsed = end - start;

   std::cout << "Matrix cost computed: " << std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count() << " ms \n";

   real_t emd_cost = 1;
   if (true) {
      matrix_t cost = compute_distance(n, 2, 1);
      if (false) {
         for (int i = 0; i < n * n; ++i) {
            for (int j = 0; j < n * n; ++j)
               fprintf(stdout, "%.2f ", cost[i][j]);
            fprintf(stdout, "\n");
         }
         exit(0);
      }

      emd_cost = compute_ns_EMD(image1, image2, cost);

      end = std::chrono::system_clock::now();
      inlineTimeElapsed = end - start;
      std::cout << "Cost: " << emd_cost << " - Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count() << " ms \n";

      start = std::chrono::system_clock::now();
   }

   auto l21_cost = apx_L2_1(image1, image2);

   // End time.
   end = std::chrono::system_clock::now();
   inlineTimeElapsed = end - start;
   std::cout << "Cost: " << l21_cost << " - Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count() << " ms \n";

   std::cout << "Ratio: " << fabs(l21_cost - emd_cost) / emd_cost * 100 << " %\n";

   return 0;
}
#endif