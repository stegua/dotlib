/**
* @fileoverview Copyright (c) 2017, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "OT_Solvers.h"

namespace DOTLib {

/**
* @brief Compute Earth Moving Distance (EMD) with L_infinity norm (with no power)
*/
real_t solve_Linf_1(const histogram_t& h1, const histogram_t& h2) {
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
   arcs.reserve(8 * d);  // at most 8 arcs for each nodes
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

   for (size_t i = 0; i < s - 1; ++i)
      for (size_t j = 0; j < s - 1; ++j) {
         arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 1, j + 1)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 1, j + 1)], nodes[ID(i, j)]));
         arcs.emplace_back(g.addArc(nodes[ID(i, j + 1)], nodes[ID(i + 1, j)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 1, j)], nodes[ID(i, j + 1)]));
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

};