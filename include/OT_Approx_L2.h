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
* @brief Compute Earth Moving Distance (EMD) between pair of images
*/
LimitValueType compute_L1_APX_EMD(const histogram_t& h1, const histogram_t& h2, const matrix_t& cost) {
   real_t distance = std::numeric_limits<real_t>::max();

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
   arcs.reserve(4 * d);
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

   fprintf(stdout, "Input graph created with %d nodes and %d arcs\n",
           countNodes(g), countArcs(g));

   NetworkSimplex<Graph, LimitValueType, LimitValueType> cycle(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g), c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i)
      b_i[nodes[i]] = int(h1[i] - h2[i]);

   // Add all edges
   for ( const auto& a : arcs ) {
      l_i[a] = 0;
      u_i[a] = cycle.INF;
      c_i[a] = 1;  // cost
   }

   //set lower/upper bounds, cost
   cycle.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   NetworkSimplex<Graph, LimitValueType, LimitValueType>::ProblemType ret = cycle.run();

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

   LimitValueType sol_value = cycle.totalCost();

   return sol_value;
}