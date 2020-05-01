/*
*  Main authors:
*     @author Stefano Gualandi <stefano.gualandi@gmail.com>
*
*     @fileoverview Copyright (c) 2017-19, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
*  Last update: April, 2019
*/

#include "DOT_BipartiteLemon.hpp"

#include <chrono>
#include <cinttypes>


// LEMON GRAPH LIBRARY
#include <lemon/smart_graph.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>

#include <lemon/cost_scaling.h>
#include <lemon/cycle_canceling.h>
#include <lemon/preflow.h>

using namespace lemon;
typedef lemon::ListDigraph LemonGraph;
typedef int64_t LimitValueType;




// Compute Kantorovich-Wasserstein distance between two measures defined in R2 (bipartite version
void DOT::BipartiteNetworkSimplex(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg) {
   int m = (int)Mu.size();
   int n = (int)Nu.size();

   // Build the graph for max flow
   LemonGraph g;

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<LemonGraph::Node> V;  // left partition
   std::vector<LemonGraph::Node> W;  // right partition

   // add first d source nodes
   for (size_t i = 0; i < m; ++i)
      V.emplace_back(g.addNode());
   for (size_t i = 0; i < n; ++i)
      W.emplace_back(g.addNode());

   std::vector<LemonGraph::Arc> arcs;
   std::vector<int64_t> a_costs;
   arcs.reserve(n*m);
   a_costs.reserve(n*m);

   for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j) {
         arcs.emplace_back(g.addArc(V[i], W[j]));
         if (algo == 0)
            a_costs.emplace_back(static_cast<LimitValueType>(DISTANCE_Rinf(Mu.getP(i), Nu.getP(j))));
         if (algo == 1)
            a_costs.emplace_back(static_cast<LimitValueType>(DISTANCE_R1(Mu.getP(i), Nu.getP(j))));
         if (algo == 2)
            a_costs.emplace_back(static_cast<LimitValueType>(DISTANCE_R2(Mu.getP(i), Nu.getP(j))));

      }

   NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<LimitValueType> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (int i = 0; i < m; ++i)
      b_i[V[i]] = +LimitValueType(Mu.getW(i));
   for (int j = 0; j < n; ++j)
      b_i[W[j]] = -LimitValueType(Nu.getW(j));

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = simplex.INF;
      c_i[a] = a_costs[i];
   }

   //set lower/upper bounds, cost
   simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   // Timinig output
   auto start = std::chrono::steady_clock::now();

   // Solve the problem to compute the distance
   NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType status = simplex.run();

   auto end = std::chrono::steady_clock::now();
   double elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   int64_t fobj = simplex.totalCost();
   fprintf(stdout, "%s %d %d Runtime %.6f Value %" PRId64 " status %d algo%d\n",
           msg.c_str(), n, countArcs(g), elapsed, fobj, status, algo);

   fflush(stdout);
}

// Compute Kantorovich-Wasserstein distance between two measures defined in R2 (bipartite version
void DOT::BipartiteCostScaling(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg) {
   int m = (int)Mu.size();
   int n = (int)Nu.size();

   // Build the graph for max flow
   LemonGraph g;

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<LemonGraph::Node> V;  // left partition
   std::vector<LemonGraph::Node> W;  // right partition

   // add first d source nodes
   for (size_t i = 0; i < m; ++i)
      V.emplace_back(g.addNode());
   for (size_t i = 0; i < n; ++i)
      W.emplace_back(g.addNode());

   std::vector<LemonGraph::Arc> arcs;
   std::vector<int64_t> a_costs;
   arcs.reserve(n*m);
   a_costs.reserve(n*m);

   for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j) {
         arcs.emplace_back(g.addArc(V[i], W[j]));
         a_costs.emplace_back(static_cast<LimitValueType>(DISTANCE_R2(Mu.getP(i), Nu.getP(j))));
      }

   CostScaling<LemonGraph, LimitValueType, LimitValueType> simplex(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<LimitValueType> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (int i = 0; i < m; ++i)
      b_i[V[i]] = +LimitValueType(Mu.getW(i));
   for (int j = 0; j < n; ++j)
      b_i[W[j]] = -LimitValueType(Nu.getW(j));

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = simplex.INF;
      c_i[a] = a_costs[i];
   }

   //set lower/upper bounds, cost
   simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   // Timinig output
   auto start = std::chrono::steady_clock::now();

   // Solve the problem to compute the distance
   CostScaling<LemonGraph, LimitValueType, LimitValueType>::ProblemType status = simplex.run();

   auto end = std::chrono::steady_clock::now();
   double elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   int64_t fobj = simplex.totalCost();
   fprintf(stdout, "%s %d %d Runtime %.6f Value %" PRId64 " status %d\n",
           msg.c_str(), n, countArcs(g), elapsed, fobj, status);

   fflush(stdout);
}

// Compute Kantorovich-Wasserstein distance between two measures defined in R2 (bipartite version
void DOT::BipartiteCycleCanceling(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg) {
   int m = (int)Mu.size();
   int n = (int)Nu.size();

   // Build the graph for max flow
   LemonGraph g;

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<LemonGraph::Node> V;  // left partition
   std::vector<LemonGraph::Node> W;  // right partition

   // add first d source nodes
   for (size_t i = 0; i < m; ++i)
      V.emplace_back(g.addNode());
   for (size_t i = 0; i < n; ++i)
      W.emplace_back(g.addNode());

   std::vector<LemonGraph::Arc> arcs;
   std::vector<LimitValueType> a_costs;
   arcs.reserve(n*m);
   a_costs.reserve(n*m);

   for (int i = 0; i < m; ++i)
      for (int j = 0; j < n; ++j) {
         arcs.emplace_back(g.addArc(V[i], W[j]));
         a_costs.emplace_back(static_cast<LimitValueType>(DISTANCE_R2(Mu.getP(i), Nu.getP(j))));
      }


   CycleCanceling<LemonGraph, LimitValueType, LimitValueType> cycle(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<LimitValueType> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (int i = 0; i < m; ++i)
      b_i[V[i]] = +LimitValueType(Mu.getW(i));
   for (int j = 0; j < n; ++j)
      b_i[W[j]] = -LimitValueType(Nu.getW(j));

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = cycle.INF;
      c_i[a] = a_costs[i];
   }

   //set lower/upper bounds, cost
   cycle.lowerMap(l_i).upperMap(u_i).costMap(c_i);
   cycle.supplyMap(b_i);

   // Timinig output
   auto start = std::chrono::steady_clock::now();

   // Solve the problem to compute the distance
   auto status = cycle.run();

   auto end = std::chrono::steady_clock::now();
   double elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   int64_t fobj = cycle.totalCost();

   fprintf(stdout, "%s %d %d Runtime %.6f Value %" PRId64 " status %d\n",
           msg.c_str(), n, countArcs(g), elapsed, fobj, status);

   fflush(stdout);
}

