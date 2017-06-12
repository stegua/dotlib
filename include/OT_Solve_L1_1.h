/**
* @fileoverview Copyright (c) 2017, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "OT_BasicTypes.h"

/**
* @brief Compute Earth Moving Distance (EMD) with L1 norm (with no power)
*        between pair of images spacing razs at "constant" angle distance, with given degree ration
*        Example: with degree=45, it puts 8 rays for vertex (at most)
*/
real_t solve_L1_1(const histogram_t& h1, const histogram_t& h2) {
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
   std::vector<real_t> arcs_costs;
   arcs.reserve(4 * d);
   for (size_t i = 0; i < s; ++i)
      for (size_t j = 0; j < s - 1; ++j) {
         arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i, j + 1)]));
         arcs.emplace_back(g.addArc(nodes[ID(i, j + 1)], nodes[ID(i, j)]));

         arcs_costs.emplace_back(1.0);
         arcs_costs.emplace_back(1.0);
      }

   for (size_t i = 0; i < s - 1; ++i)
      for (size_t j = 0; j < s; ++j) {
         arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 1, j)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 1, j)], nodes[ID(i, j)]));

         arcs_costs.emplace_back(1.0);
         arcs_costs.emplace_back(1.0);
      }

   for (size_t i = 0; i < s - 1; ++i)
      for (size_t j = 0; j < s - 1; ++j) {
         arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 1, j + 1)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 1, j + 1)], nodes[ID(i, j)]));
         arcs.emplace_back(g.addArc(nodes[ID(i, j + 1)], nodes[ID(i + 1, j)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 1, j)], nodes[ID(i, j + 1)]));

         arcs_costs.emplace_back(sqrt(2.0));
         arcs_costs.emplace_back(sqrt(2.0));
         arcs_costs.emplace_back(sqrt(2.0));
         arcs_costs.emplace_back(sqrt(2.0));
      }

   fprintf(stdout, "Input graph created with %d nodes and %d arcs\n",
           countNodes(g), countArcs(g));

   NetworkSimplex<Graph, LimitValueType, real_t> cycle(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<real_t> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i)
      b_i[nodes[i]] = int(h1[i] - h2[i]);

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = cycle.INF;
      c_i[a] = std::trunc(65536 * arcs_costs[i]) / 65536;
   }

   //set lower/upper bounds, cost
   cycle.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

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