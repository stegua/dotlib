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

void flowCycleCanceling() {
   /**
   * calculate min cost max flow via cycle canceling
   * sample graph (arcs are dircted from left to right
   *
   *     2     5
   *    / \   / \
   *   /   \ /   \
   *  0     4--6--1
   *   \   / \   /
   *    \ /   \ /
   *     3     7
   *
   * vertex 0 represents the source, 1 target
   * cost are set to 1, ecept for source and target arcs
   * cost for (3,4) = 2
   */

   Graph g;

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g), c_i(g);

   // source and sink vertices
   Graph::Node n0 = g.addNode(), n1 = g.addNode();

   // other nodes/vertices
   Graph::Node n2 = g.addNode(), n3 = g.addNode(), n4 = g.addNode(), n5 = g.addNode(),
               n6 = g.addNode(), n7 = g.addNode();

   Graph::Arc a02 = g.addArc(n0, n2), a03 = g.addArc(n0, n3), a24 = g.addArc(n2, n4),
              a34 = g.addArc(n3, n4), a45 = g.addArc(n4, n5), a46 = g.addArc(n4, n6),
              a47 = g.addArc(n4, n7), a51 = g.addArc(n5, n1), a61 = g.addArc(n6, n1),
              a71 = g.addArc(n7, n1);

   std::vector<Graph::Arc> arcs = { a02, a03, a24, a34, a45, a46, a47, a51, a61, a71 };

   // arcs for circulation! (source to target and vice versa)
   Graph::Arc ts = g.addArc(n1, n0);
   Graph::Arc st = g.addArc(n0, n1);


   // graph must be finished before initializing cycle canceling !
   CycleCanceling<Graph, LimitValueType, LimitValueType> cycle(g);

   // bounds for circulation arcs
   l_i[ts] = 0;
   u_i[ts] = cycle.INF;
   c_i[ts] = 0;
   l_i[st] = 0;
   u_i[ts] = cycle.INF;
   c_i[ts] = 0;

   // set cost for arcs
   for (auto & a : arcs) {
      l_i[a] = 1;
      u_i[a] = cycle.INF;
      c_i[a] = (g.source(a) == n0 || g.target(a) == n1) ? 0 : 1;
   }
   c_i[a34] = 2; // arc (3,4)

   //set lower/upper bounds, cost
   cycle.lowerMap(l_i).upperMap(u_i).costMap(c_i);

   CycleCanceling<Graph, LimitValueType>::ProblemType ret = cycle.run();

   // get flow of arcs
   ListDigraph::ArcMap<LimitValueType> flow(g);
   cycle.flowMap(flow);

   digraphWriter(g).                  // write g to the standard output
   arcMap("cost", c_i).          // write 'cost' for for arcs
   arcMap("flow", flow).          // write 'flow' for for arcs
   arcMap("l_i", l_i).
   arcMap("u_i", u_i).
   node("source", n0).             // write s to 'source'
   node("target", n1).             // write t to 'target'
   run();

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

   // determine flow from source and flow to target
   double flowToT = 0, flowToS = 0;
   for (ListDigraph::ArcMap<LimitValueType>::ItemIt a(flow); a != INVALID; ++a) {
      if (g.target(a) == n1)
         flowToT += cycle.flow(a);

      if (g.source(a) == n0)
         flowToS += cycle.flow(a);
   }

   std::cerr << "flow to t = " << flowToT << " flow from s = " << flowToS << std::endl;
   std::cerr << "total cost = " << cycle.totalCost<double>() << std::endl;

}

/**
* @brief Read images from text file and return corresponfing matrix
*/
histogram_t read_image(size_t n, std::string filename) {
   auto m = [&n](size_t x, size_t y) {
      return x * n + y;
   };

   histogram_t image(n * n, 0);
   image.shrink_to_fit();

   std::ifstream in_file(filename);

   if (!in_file) {
      fprintf(stdout, "Cannot open file %s.\n", filename);
      return image;
   }

   for (size_t i = 0; i < n; ++i) {
      size_t j = 0;
      std::string         line;
      std::getline(in_file, line);
      std::stringstream   lineStream(line);
      std::string         cell;

      while (std::getline(lineStream, cell, ',')) {
         image[m(i, j)] = stoi(cell);
         ++j;
      }
   }

   in_file.close();

   return image;
}

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
               cost[m(i, j)][m(v, w)] = pow(pow(pow(fabs(i - v), p) + pow(fabs(j - w), p), 1.0 / p), q);
      }

   return cost;
}

/**
* @brief Compute Earth Moving Distance (EMD) between pair of images
*/
int compute_EMD(const histogram_t& h1, const histogram_t& h2, const matrix_t& cost) {
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

int main(int argc, char* argv) {
   // Size of the square image
   size_t n = 32;

   histogram_t image1 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data32_1001.csv");
   histogram_t image2 = read_image(n, "D:\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data32_1010.csv");

   if (false) {
      n = 2;
      image1 = { 1, 0, 0, 0 };
      image2 = { 0, 0, 0, 1 };
   }

   // Time vars
   std::chrono::time_point<std::chrono::system_clock> start, end;
   // Start time.

   start = std::chrono::system_clock::now();

   matrix_t cost = compute_distance(n, 2, 1);
   if (false) {
      for (int i = 0; i < n * n; ++i) {
         for (int j = 0; j < n * n; ++j)
            fprintf(stdout, "%.2f ", cost[i][j]);
         fprintf(stdout, "\n");
      }
      exit(0);
   }

   // End time.
   end = std::chrono::system_clock::now();
   std::chrono::duration<double> inlineTimeElapsed = end - start;

   std::cout << "Matrix cost copmuted: " << std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count() << " ms \n";

   real_t emd_cost = compute_ns_EMD(image1, image2, cost);

   end = std::chrono::system_clock::now();
   inlineTimeElapsed = end - start;
   std::cout << "Cost: " << emd_cost << " - Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count() << " ms \n";

   start = std::chrono::system_clock::now();
   auto l21_cost = apx_L2_1(image1, image2);
   std::cout << "Cost: " << l21_cost << " - Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(inlineTimeElapsed).count() << " ms \n";

   std::cout << "Ratio: " << fabs(l21_cost - emd_cost) / emd_cost * 100 << " %\n";

   // End time.
   end = std::chrono::system_clock::now();
   inlineTimeElapsed = end - start;


   return 0;
}