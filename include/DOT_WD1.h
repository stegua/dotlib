/**
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "DOT_BasicTypes.h"
#include "DOT_Histogram2D.h"

namespace DOT {

/**
* @brief Compute the Wasserstein distance of order 1 with ground distance L1
*/
int64_t solve_network_L1(const Histogram2D& h1, const Histogram2D& h2) {
   // Build the graph for max flow
   LemonGraph g;

   size_t s = h1.getN();
   size_t d = s * s;

   auto ID = [&s](size_t x, size_t y) {
      return x * s + y;
   };

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<LemonGraph::Node> nodes;
   // add first d source nodes
   for (size_t i = 0; i < d; ++i)
      nodes.emplace_back(g.addNode());

   std::vector<LemonGraph::Arc> arcs;
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

   fprintf(stdout, "Input graph created with %ld nodes and %ld arcs\n", countNodes(g), countArcs(g));

   LemonSimplex simplex(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<LimitValueType> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   {
      size_t idx = 0;
      for (size_t i = 0; i < s; ++i)
         for (size_t j = 0; j < s; ++j)
            b_i[nodes[idx++]] = LimitValueType(h1.get(i, j) - h2.get(i, j));
   }

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = simplex.INF;
      c_i[a] = 1;
   }

   //set lower/upper bounds, cost
   simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   // Solve the problem to compute the distance
   NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType ret = simplex.run();

   switch (ret) {
   case NetworkSimplex<LemonGraph>::INFEASIBLE:
      fprintf(stdout, "NetworkSimplex<Graph>::INFEASIBLE\n");
      break;
   case NetworkSimplex<LemonGraph>::OPTIMAL:
      fprintf(stdout, "NetworkSimplex<Graph>::OPTIMAL\n");
      break;
   case NetworkSimplex<LemonGraph>::UNBOUNDED:
      fprintf(stdout, "NetworkSimplex<Graph>::UNBOUNDED\n");
      break;
   }

   return simplex.totalCost();
}

}