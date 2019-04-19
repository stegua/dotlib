/*
*  Main authors:
*     @author Stefano Gualandi <stefano.gualandi@gmail.com>
*
*     @fileoverview Copyright (c) 2017-19, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
*  Last update: April, 2019
*/

#include "DOT_KWD_P1.hpp"

#include <chrono>
#include <cassert>
#include <cinttypes>

// LEMON GRAPH LIBRARY
#include <lemon/smart_graph.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>

#include <lemon/cost_scaling.h>

using namespace lemon;
typedef lemon::ListDigraph LemonGraph;
typedef int64_t LimitValueType;


// Compute EXACT Kantorovich-Wasserstein Distance (KWD) of order 1,
// using the BGV method with ground distance L1 (Manhattan)
void DOT::BGV_L1(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg) {
   int m = (int)Mu.size();
   int n = (int)Nu.size();

   assert(n == m);

   size_t d = m;
   size_t s = static_cast<size_t>(sqrt(d));

   auto ID = [&s](size_t x, size_t y) {
      return x * s + y;
   };

   // Build the graph for max flow
   LemonGraph g;

   // add d nodes
   std::vector<LemonGraph::Node> V;  // left partition

   // add first d source nodes
   for (size_t i = 0; i < d; ++i)
      V.emplace_back(g.addNode());

   std::vector<LemonGraph::Arc> arcs;
   arcs.reserve(4 * d); // at most 4 arcs for each nodes

   for (size_t i = 0; i < s; ++i)
      for (size_t j = 0; j < s - 1; ++j) {
         arcs.emplace_back(g.addArc(V[ID(i, j)], V[ID(i, j + 1)]));
         arcs.emplace_back(g.addArc(V[ID(i, j + 1)], V[ID(i, j)]));
      }

   for (size_t i = 0; i < s - 1; ++i)
      for (size_t j = 0; j < s; ++j) {
         arcs.emplace_back(g.addArc(V[ID(i, j)], V[ID(i + 1, j)]));
         arcs.emplace_back(g.addArc(V[ID(i + 1, j)], V[ID(i, j)]));
      }

   arcs.shrink_to_fit();

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<LimitValueType> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (int i = 0; i < d; ++i)
      b_i[V[i]] = LimitValueType(Mu.getW(i) - Nu.getW(i));

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = std::numeric_limits<LimitValueType>::max();
      c_i[a] = 1;
   }

   // Timinig output
   auto start = std::chrono::steady_clock::now();
   int64_t fobj = 0;
   int status = 0;

   if (algo == 0) {
      //set lower/upper bounds, cost
      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);
      simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

      // Solve the problem to compute the distance
      status = simplex.run();

      fobj = simplex.totalCost();
   } else {
      //set lower/upper bounds, cost
      CostScaling<LemonGraph, LimitValueType, LimitValueType> scaling(g);
      scaling.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

      // Solve the problem to compute the distance
      status = scaling.run();

      fobj = scaling.totalCost();
   }

   auto end = std::chrono::steady_clock::now();
   double elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   fprintf(stdout, "%s %d %d Runtime %.6f Value %" PRId64 " status %d\n",
           msg.c_str(), n, countArcs(g), elapsed, fobj, status);

   fflush(stdout);
}

// Compute EXACT Kantorovich-Wasserstein Distance (KWD) of order 1,
// using the BGV method with ground distance L_inf
void DOT::BGV_Linf(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg) {
   int m = (int)Mu.size();
   int n = (int)Nu.size();

   assert(n == m);

   size_t d = m;
   size_t s = static_cast<size_t>(sqrt(d));

   auto ID = [&s](size_t x, size_t y) {
      return x * s + y;
   };

   // Build the graph for max flow
   LemonGraph g;

   // add d nodes
   std::vector<LemonGraph::Node> V;  // left partition

   // add first d source nodes
   for (size_t i = 0; i < d; ++i)
      V.emplace_back(g.addNode());

   std::vector<LemonGraph::Arc> arcs;
   arcs.reserve(8 * d);

   for (size_t i = 0; i < s; ++i)
      for (size_t j = 0; j < s - 1; ++j) {
         arcs.emplace_back(g.addArc(V[ID(i, j)], V[ID(i, j + 1)]));
         arcs.emplace_back(g.addArc(V[ID(i, j + 1)], V[ID(i, j)]));
      }

   for (size_t i = 0; i < s - 1; ++i)
      for (size_t j = 0; j < s; ++j) {
         arcs.emplace_back(g.addArc(V[ID(i, j)], V[ID(i + 1, j)]));
         arcs.emplace_back(g.addArc(V[ID(i + 1, j)], V[ID(i, j)]));
      }

   for (size_t i = 0; i < s - 1; ++i)
      for (size_t j = 0; j < s - 1; ++j) {
         arcs.emplace_back(g.addArc(V[ID(i, j)], V[ID(i + 1, j + 1)]));
         arcs.emplace_back(g.addArc(V[ID(i + 1, j + 1)], V[ID(i, j)]));
         arcs.emplace_back(g.addArc(V[ID(i, j + 1)], V[ID(i + 1, j)]));
         arcs.emplace_back(g.addArc(V[ID(i + 1, j)], V[ID(i, j + 1)]));
      }

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<LimitValueType> c_i(g);

   // Flow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i)
      b_i[V[i]] = LimitValueType(Mu.getW(i) - Nu.getW(i));

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = std::numeric_limits<LimitValueType>::max();
      c_i[a] = 1;
   }

   // Timinig output
   auto start = std::chrono::steady_clock::now();
   int64_t fobj = 0;
   int status = 0;

   if (algo == 0) {
      //set lower/upper bounds, cost
      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);
      simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

      // Solve the problem to compute the distance
      status = simplex.run();

      fobj = simplex.totalCost();
   } else {
      //set lower/upper bounds, cost
      CostScaling<LemonGraph, LimitValueType, LimitValueType> scaling(g);
      scaling.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

      // Solve the problem to compute the distance
      status = scaling.run();

      fobj = scaling.totalCost();
   }

   auto end = std::chrono::steady_clock::now();
   double elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   fprintf(stdout, "%s %d %d Runtime %.6f Value %" PRId64 " status %d\n",
           msg.c_str(), n, countArcs(g), elapsed, fobj, status);

   fflush(stdout);
}

// Compute EXACT Kantorovich-Wasserstein Distance (KWD) of order 1,
// using the BGV method with ground distance L2
void DOT::BGV_L2(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg,
                 const std::vector<std::pair<int, int>>& coprimes) {
   int m = (int)Mu.size();
   int n = (int)Nu.size();

   assert(n == m);

   size_t d = m;
   int s = static_cast<int>(sqrt(double(d)));

   auto ID = [&s](size_t x, size_t y) {
      return x * s + y;
   };

   // Build the graph for max flow
   LemonGraph g;

   // add d nodes
   std::vector<LemonGraph::Node> V;  // left partition

   // add first d source nodes
   for (size_t i = 0; i < d; ++i)
      V.emplace_back(g.addNode());

   std::vector<LemonGraph::Arc> arcs;

   std::vector<LimitValueType> a_costs;

   for (int i = 0; i < s; ++i)
      for (int j = 0; j < s; ++j)
         for (const auto& p : coprimes) {
            int v = p.first;
            int w = p.second;
            if (i + v >= 0 && i + v < s && j + w >= 0 && j + w < s) {
               arcs.emplace_back(g.addArc(V[ID(i, j)], V[ID(i + v, j + w)]));
               a_costs.emplace_back(static_cast<LimitValueType>(DISTANCE_R2(Mu.getP(ID(i, j)), Mu.getP(ID(i + v, j + w)))));
            }
         }

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<LimitValueType> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (int i = 0; i < d; ++i)
      b_i[V[i]] = LimitValueType(Mu.getW(i) - Nu.getW(i));

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = std::numeric_limits<LimitValueType>::max();
      c_i[a] = a_costs[i];
   }

   NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);

   //set lower/upper bounds, cost
   simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   // Timinig output
   auto start = std::chrono::steady_clock::now();

   // Solve the problem to compute the distance
   NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType status = simplex.run();

   auto end = std::chrono::steady_clock::now();
   double elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   int64_t fobj = simplex.totalCost();

   fprintf(stdout, "%s %d %d Runtime %.6f Value %" PRId64 " status %d\n",
           msg.c_str(), n, countArcs(g), elapsed, fobj, status);

   fflush(stdout);
}