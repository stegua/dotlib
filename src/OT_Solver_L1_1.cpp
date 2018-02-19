/**
* @fileoverview Copyright (c) 2017, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include <numeric>
#include "OT_Solvers.h"

namespace DOTLib {

/**
* @brief Compute Earth Moving Distance (EMD) with L1 norm (with no power)
*        between pair of images spacing razs at "constant" angle distance, with given degree ration
*        Example: with degree=45, it puts 8 rays for vertex (at most)
*/
real_t solve_L1_1(const histogram_t& h1, const histogram_t& h2) {
   auto logger = spd::get("console");

   size_t d = h1.size();
   size_t s = static_cast<size_t>(sqrt(d));

   auto ID = [&s](size_t x, size_t y) {
      return x * s + y;
   };

   // Build the graph for max flow
   Graph g;

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<Graph::Node> nodes;
   // add first d source nodes
   for (size_t i = 0; i < d; ++i)
      nodes.emplace_back(g.addNode());

   std::vector<Graph::Arc> arcs;
   arcs.reserve(4 * d); // at most 4 arcs for each nodes
   for (size_t i = 0; i < s; ++i)
      for (size_t j = 0; j < s - 1; ++j) {
         arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i, j + 1)]));
         arcs.emplace_back(g.addArc(nodes[ID(i, j + 1)], nodes[ID(i, j)]));
      }

   for (size_t i = 0; i < s - 1; ++i)
      for (size_t j = 0; j < s; ++j) {
         arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 1, j)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 1, j)], nodes[ID(i, j)]));
      }

   logger->info("Input graph created with {1} nodes and {0} arcs", countArcs(g), countNodes(g));

   NetworkSimplex<Graph, LimitValueType, real_t> simplex(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<real_t> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i)
      b_i[nodes[i]] = LimitValueType(h1[i] - h2[i]);

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = simplex.INF;
      c_i[a] = 1.0;
   }

   //set lower/upper bounds, cost
   simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   NetworkSimplex<Graph, LimitValueType, real_t>::ProblemType ret = simplex.run();

   switch (ret) {
   case NetworkSimplex<Graph>::INFEASIBLE:
      logger->error("NetworkSimplex<Graph>::INFEASIBLE");
      break;
   case NetworkSimplex<Graph>::OPTIMAL:
      logger->info("NetworkSimplex<Graph>::OPTIMAL");
      break;
   case NetworkSimplex<Graph>::UNBOUNDED:
      logger->error("NetworkSimplex<Graph>::UNBOUNDED");
      break;
   }

   real_t sol_value = simplex.totalCost();

   return sol_value;
}


/**
* @brief Compute Earth Moving Distance (EMD) with L1 norm (with no power)
*        between pair of images spacing razs at "constant" angle distance, with given degree ration
*        Example: with degree=45, it puts 8 rays for vertex (at most)
*/
real_t solve_L2_1(const histogram_t& h1, const histogram_t& h2, int L) {
   auto logger = spd::get("logs");

   size_t d = h1.size();
   size_t s = static_cast<size_t>(sqrt(d));

   auto ID = [&s](size_t x, size_t y) {
      return x * s + y;
   };

   // Time vars
   std::chrono::time_point<std::chrono::system_clock> start, end;
   // Start time.
   start = std::chrono::system_clock::now();

   // Build the graph for max flow
   Graph g;
   g.reserveArc(577 * d);
   g.reserveNode(d);

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<Graph::Node> nodes;
   nodes.reserve(d);
   // add first d source nodes
   for (size_t i = 0; i < d; ++i)
      nodes.emplace_back(g.addNode());

   std::vector<Graph::Arc> arcs;
   std::vector<real_t> a_costs;
   arcs.reserve(577 *d); // at most 4 arcs for each nodes
   a_costs.reserve(577 * d);

   std::vector<std::pair<int, int>> AL;
   for (int v = -L; v <= L; ++v)
      for (int w = -L; w <= L; ++w)
         if (!(v == 0 && w == 0) && std::gcd(v, w) == 1) {
            AL.emplace_back(v, w);
         }

   for (int i = 0; i < s; ++i)
      for (int j = 0; j < s; ++j)
         for (const auto& p: AL) {
            int v = p.first;
            int w = p.second;
            if (i + v >= 0 && i + v < s && j + w >= 0 && j + w < s) {
               arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + v, j + w)]));
               a_costs.emplace_back(sqrt(pow(v, 2) + pow(w, 2)));
            }
         }

   // End time.
   end = std::chrono::system_clock::now();
   std::chrono::duration<double> inlineTimeElapsed = end - start;

   start = std::chrono::system_clock::now();

   NetworkSimplex<Graph, LimitValueType, real_t> simplex(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i)
      b_i[nodes[i]] = LimitValueType(h1[i] - h2[i]);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<real_t> c_i(g);

   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = simplex.INF;
      c_i[a] = std::trunc(MUL*a_costs[i]);
   }

   //set lower/upper bounds, cost
   simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   NetworkSimplex<Graph, LimitValueType, real_t>::ProblemType ret = simplex.run();

   switch (ret) {
   case NetworkSimplex<Graph>::INFEASIBLE:
      logger->error("NetworkSimplex<Graph>::INFEASIBLE");
      break;
   case NetworkSimplex<Graph>::OPTIMAL:
      logger->info("NetworkSimplex<Graph>::OPTIMAL");
      break;
   case NetworkSimplex<Graph>::UNBOUNDED:
      logger->error("NetworkSimplex<Graph>::UNBOUNDED");
      break;
   }

   real_t sol_value = simplex.totalCost()/ MUL;

   end = std::chrono::system_clock::now();
   std::chrono::duration<double> run_time = end - start;

   logger->info("STE nodes {} arcs {} n {} N {} build_time {} cost {:.2f} run_time {} L {}",
                countNodes(g), countArcs(g), h1.size(), s, inlineTimeElapsed.count(), sol_value, run_time.count(), L);

   return sol_value;
}

};