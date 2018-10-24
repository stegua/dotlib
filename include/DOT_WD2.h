/**
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "DOT_BasicTypes.h"
#include "DOT_Histogram.h"
#include "DOT_Histogram2D.h"

// LEMON GRAPH LIBRARY
#include <lemon/smart_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <lemon/cost_scaling.h>
#include <lemon/cycle_canceling.h>

#include <lemon/preflow.h>

using namespace lemon;
typedef lemon::ListDigraph LemonGraph;
typedef int64_t LimitValueType;
typedef NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> LemonSimplex;


namespace DOT {

/**
* @brief Compute the Wasserstein distance of order 2 with ground distance L2 (Euclidean norm)
*/
int64_t solve_bipartite_network_L2(const Histogram2D& h1, const Histogram2D& h2) {
   // Build the graph for max flow
   LemonGraph g;

   int s = static_cast<int>(h1.getN());
   size_t d = s * s;

   auto ID = [&s](size_t x, size_t y) {
      return x * s + y;
   };

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<LemonGraph::Node> V;  // left partition
   std::vector<LemonGraph::Node> W;  // right partition

   // add first d source nodes
   for (size_t i = 0; i < d; ++i)
      V.emplace_back(g.addNode());
   for (size_t i = 0; i < d; ++i)
      W.emplace_back(g.addNode());

   std::vector<LemonGraph::Arc> arcs;
   std::vector<int64_t> a_costs;

   for (int i = 0; i < s; ++i)
      for (int j = 0; j < s; ++j)
         for (int v = 0; v < s; ++v)
            for (int w = 0; w < s; ++w) {
               arcs.emplace_back(g.addArc(V[ID(i, j)], W[ID(v, w)]));
               a_costs.emplace_back(static_cast<int64_t>(pow(i - v, 2) + pow(j - w, 2)));
            }

   logger.info("Input graph created with %d nodes and %d arcs", countNodes(g), countArcs(g));

   NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<LimitValueType> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   {
      size_t idx = 0;
      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j) {
            b_i[V[ID(i, j)]] = +LimitValueType(h1.get(i, j));
            b_i[W[ID(i, j)]] = -LimitValueType(h2.get(i, j));
         }
   }

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = simplex.INF;
      c_i[a] = a_costs[i];
   }

   //set lower/upper bounds, cost
   simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   // Solve the problem to compute the distance
   NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType ret = simplex.run();

   switch (ret) {
   case NetworkSimplex<LemonGraph>::INFEASIBLE:
      logger.info("NetworkSimplex<Graph>::INFEASIBLE");
      break;
   case NetworkSimplex<LemonGraph>::OPTIMAL:
      logger.info("NetworkSimplex<Graph>::OPTIMAL");
      break;
   case NetworkSimplex<LemonGraph>::UNBOUNDED:
      logger.info("NetworkSimplex<Graph>::UNBOUNDED");
      break;
   }

   return simplex.totalCost();
}

/**
* @brief Compute the Wasserstein distance of order 2 with ground distance L2 (Euclidean norm)
*/
int64_t solve_bipartite_network_L2_fcs(const Histogram& h1, const Histogram& h2) {
   // Build the graph for max flow
   LemonGraph g;

   int s = static_cast<int>(h1.getN());
   if (h1.getD() == 2) {
      size_t d = s * s;

      auto ID = [&s](size_t x, size_t y) {
         return x * s + y;
      };

      // add d nodes for each histrogam (d+1) source, (d+2) target
      std::vector<LemonGraph::Node> V;  // left partition
      std::vector<LemonGraph::Node> W;  // right partition

      // add first d source nodes
      for (size_t i = 0; i < d; ++i)
         V.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         W.emplace_back(g.addNode());

      std::vector<LemonGraph::Arc> arcs;
      std::vector<int64_t> a_costs;

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int v = 0; v < s; ++v)
               for (int w = 0; w < s; ++w) {
                  arcs.emplace_back(g.addArc(V[ID(i, j)], W[ID(v, w)]));
                  a_costs.emplace_back(static_cast<int64_t>(pow(i - v, 2) + pow(j - w, 2)));
               }

      logger.info("Input graph created with %d nodes and %d arcs", countNodes(g), countArcs(g));

      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);

      // lower and upper bounds, cost
      ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
      ListDigraph::ArcMap<LimitValueType> c_i(g);

      // FLow balance
      ListDigraph::NodeMap<LimitValueType> b_i(g);
      {
         size_t idx = 0;
         for (int i = 0; i < s; ++i)
            for (int j = 0; j < s; ++j) {
               b_i[V[ID(i, j)]] = +LimitValueType(h1.get(i, j));
               b_i[W[ID(i, j)]] = -LimitValueType(h2.get(i, j));
            }
      }

      // Add all edges
      for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
         const auto& a = arcs[i];
         l_i[a] = 0;
         u_i[a] = simplex.INF;
         c_i[a] = a_costs[i];
      }

      //set lower/upper bounds, cost
      simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

      // Solve the problem to compute the distance
      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType ret = simplex.run();

      switch (ret) {
      case NetworkSimplex<LemonGraph>::INFEASIBLE:
         logger.info("NetworkSimplex<Graph>::INFEASIBLE");
         break;
      case NetworkSimplex<LemonGraph>::OPTIMAL:
         logger.info("NetworkSimplex<Graph>::OPTIMAL");
         break;
      case NetworkSimplex<LemonGraph>::UNBOUNDED:
         logger.info("NetworkSimplex<Graph>::UNBOUNDED");
         break;
      }

      return simplex.totalCost();
   }
   if (h1.getD() == 3) {
      size_t d = s * s * s;

      auto ID = [&s](size_t x, size_t y, size_t z) {
         return x * s*s + y*s + z;
      };

      // add d nodes for each histrogam (d+1) source, (d+2) target
      std::vector<LemonGraph::Node> V;  // left partition
      std::vector<LemonGraph::Node> W;  // right partition

      // add first d source nodes
      for (size_t i = 0; i < d; ++i)
         V.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         W.emplace_back(g.addNode());

      std::vector<LemonGraph::Arc> arcs;
      std::vector<int64_t> a_costs;

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int h = 0; h < s; ++h)
               for (int v = 0; v < s; ++v)
                  for (int w = 0; w < s; ++w)
                     for (int u = 0; u < s; ++u) {
                        arcs.emplace_back(g.addArc(V[ID(i, j, h)], W[ID(v, w, u)]));
                        a_costs.emplace_back(static_cast<int64_t>(pow(i - v, 2) + pow(j - w, 2) + pow(h - u, 2)));
                     }

      logger.info("Input graph created with %d nodes and %d arcs", countNodes(g), countArcs(g));

      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);

      // lower and upper bounds, cost
      ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
      ListDigraph::ArcMap<LimitValueType> c_i(g);

      // FLow balance
      ListDigraph::NodeMap<LimitValueType> b_i(g);
      {
         size_t idx = 0;
         for (int i = 0; i < s; ++i)
            for (int j = 0; j < s; ++j)
               for (int h = 0; h < s; ++h) {
                  b_i[V[ID(i, j, h)]] = +LimitValueType(h1.get(i, j, h));
                  b_i[W[ID(i, j, h)]] = -LimitValueType(h2.get(i, j, h));
               }
      }

      // Add all edges
      for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
         const auto& a = arcs[i];
         l_i[a] = 0;
         u_i[a] = simplex.INF;
         c_i[a] = a_costs[i];
      }

      //set lower/upper bounds, cost
      simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

      // Solve the problem to compute the distance
      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType ret = simplex.run();

      switch (ret) {
      case NetworkSimplex<LemonGraph>::INFEASIBLE:
         logger.info("NetworkSimplex<Graph>::INFEASIBLE");
         break;
      case NetworkSimplex<LemonGraph>::OPTIMAL:
         logger.info("NetworkSimplex<Graph>::OPTIMAL");
         break;
      case NetworkSimplex<LemonGraph>::UNBOUNDED:
         logger.info("NetworkSimplex<Graph>::UNBOUNDED");
         break;
      }

      return simplex.totalCost();
   }

   if (h1.getD() == 4) {
      size_t d = s * s * s * s;

      auto ID = [&s](size_t x, size_t y, size_t z, size_t u) {
         return x * s*s*s + y * s*s + z*s + u;
      };

      // add d nodes for each histrogam (d+1) source, (d+2) target
      std::vector<LemonGraph::Node> V;  // left partition
      std::vector<LemonGraph::Node> W;  // right partition

      // add first d source nodes
      for (size_t i = 0; i < d; ++i)
         V.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         W.emplace_back(g.addNode());

      std::vector<LemonGraph::Arc> arcs;
      std::vector<int64_t> a_costs;

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int h = 0; h < s; ++h)
               for (int l = 0; l < s; ++l)
                  for (int v = 0; v < s; ++v)
                     for (int w = 0; w < s; ++w)
                        for (int u = 0; u < s; ++u)
                           for (int z = 0; z < s; ++z) {
                              arcs.emplace_back(g.addArc(V[ID(i, j, h, l)], W[ID(v, w, u, z)]));
                              a_costs.emplace_back(static_cast<int64_t>(pow(i - v, 2) + pow(j - w, 2) + pow(h - u, 2) + pow(l - z, 2)));
                           }

      logger.info("Input graph created with %d nodes and %d arcs", countNodes(g), countArcs(g));

      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);

      // lower and upper bounds, cost
      ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
      ListDigraph::ArcMap<LimitValueType> c_i(g);

      // FLow balance
      ListDigraph::NodeMap<LimitValueType> b_i(g);
      {
         size_t idx = 0;
         for (int i = 0; i < s; ++i)
            for (int j = 0; j < s; ++j)
               for (int h = 0; h < s; ++h)
                  for (int l = 0; l < s; ++l) {
                     b_i[V[ID(i, j, h, l)]] = +LimitValueType(h1.get(i, j, h, l));
                     b_i[W[ID(i, j, h, l)]] = -LimitValueType(h2.get(i, j, h, l));
                  }
      }

      // Add all edges
      for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
         const auto& a = arcs[i];
         l_i[a] = 0;
         u_i[a] = simplex.INF;
         c_i[a] = a_costs[i];
      }

      //set lower/upper bounds, cost
      simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

      // Solve the problem to compute the distance
      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType ret = simplex.run();

      switch (ret) {
      case NetworkSimplex<LemonGraph>::INFEASIBLE:
         logger.info("NetworkSimplex<Graph>::INFEASIBLE");
         break;
      case NetworkSimplex<LemonGraph>::OPTIMAL:
         logger.info("NetworkSimplex<Graph>::OPTIMAL");
         break;
      case NetworkSimplex<LemonGraph>::UNBOUNDED:
         logger.info("NetworkSimplex<Graph>::UNBOUNDED");
         break;
      }

      return simplex.totalCost();
   }

}

/**
* @brief Compute the Wasserstein distance of order 2 with ground distance L2 (Euclidean norm)
*/
int64_t solve_tripartite_network_L2(const Histogram2D& h1, const Histogram2D& h2) {
   // Build the graph for max flow
   LemonGraph g;

   int s = static_cast<int>(h1.getN());

   size_t d = s * s;

   auto ID = [&s](size_t x, size_t y) {
      return x * s + y;
   };

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<LemonGraph::Node> U;  // left partition
   std::vector<LemonGraph::Node> V;  // middle partition
   std::vector<LemonGraph::Node> W;  // right partition

   // add first d source nodes
   for (size_t i = 0; i < d; ++i)
      U.emplace_back(g.addNode());
   for (size_t i = 0; i < d; ++i)
      V.emplace_back(g.addNode());
   for (size_t i = 0; i < d; ++i)
      W.emplace_back(g.addNode());

   std::vector<LemonGraph::Arc> arcs;
   std::vector<int64_t> a_costs;

   for (int i = 0; i < s; ++i)
      for (int j = 0; j < s; ++j)
         for (int v = 0; v < s; ++v) {
            arcs.emplace_back(g.addArc(U[ID(i, j)], V[ID(i, v)]));
            a_costs.emplace_back(static_cast<int64_t>(pow(j - v, 2)));
         }

   for (int i = 0; i < s; ++i)
      for (int j = 0; j < s; ++j)
         for (int v = 0; v < s; ++v) {
            arcs.emplace_back(g.addArc(V[ID(i, j)], W[ID(v, j)]));
            a_costs.emplace_back(static_cast<int64_t>(pow(i - v, 2)));
         }

   logger.info("Input graph created with %d nodes and %d arcs", countNodes(g), countArcs(g));

   NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<LimitValueType> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   {
      size_t idx = 0;
      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j) {
            b_i[U[ID(i, j)]] = +LimitValueType(h1.get(i, j));
            b_i[V[ID(i, j)]] =  LimitValueType(0);
            b_i[W[ID(i, j)]] = -LimitValueType(h2.get(i, j));
         }
   }

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = simplex.INF;
      c_i[a] = a_costs[i];
   }

   //set lower/upper bounds, cost
   simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   // Solve the problem to compute the distance
   NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType ret = simplex.run();

   switch (ret) {
   case NetworkSimplex<LemonGraph>::INFEASIBLE:
      logger.info("NetworkSimplex<Graph>::INFEASIBLE");
      break;
   case NetworkSimplex<LemonGraph>::OPTIMAL:
      logger.info("NetworkSimplex<Graph>::OPTIMAL");
      break;
   case NetworkSimplex<LemonGraph>::UNBOUNDED:
      logger.info("NetworkSimplex<Graph>::UNBOUNDED");
      break;
   }

   return simplex.totalCost();
}

/**
* @brief Compute the Wasserstein distance of order 2 with ground distance L2 (Euclidean norm)
*/
int64_t solve_tripartite_network_L2_fcs(const Histogram& h1, const Histogram& h2) {
   // Build the graph for max flow
   LemonGraph g;

   int s = static_cast<int>(h1.getN());
   if (h1.getD() == 2) {
      size_t d = s * s;

      auto ID = [&s](size_t x, size_t y) {
         return x * s + y;
      };

      // add d nodes for each histrogam (d+1) source, (d+2) target
      std::vector<LemonGraph::Node> U;  // left partition
      std::vector<LemonGraph::Node> V;  // middle partition
      std::vector<LemonGraph::Node> W;  // right partition

      // add first d source nodes
      for (size_t i = 0; i < d; ++i)
         U.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         V.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         W.emplace_back(g.addNode());

      std::vector<LemonGraph::Arc> arcs;
      std::vector<int64_t> a_costs;

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int v = 0; v < s; ++v) {
               arcs.emplace_back(g.addArc(U[ID(i, j)], V[ID(i, v)]));
               a_costs.emplace_back(static_cast<int64_t>(pow(j - v, 2)));
            }

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int v = 0; v < s; ++v) {
               arcs.emplace_back(g.addArc(V[ID(i, j)], W[ID(v, j)]));
               a_costs.emplace_back(static_cast<int64_t>(pow(i - v, 2)));
            }

      logger.info("Input graph created with %d nodes and %d arcs", countNodes(g), countArcs(g));

      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);

      // lower and upper bounds, cost
      ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
      ListDigraph::ArcMap<LimitValueType> c_i(g);

      // FLow balance
      ListDigraph::NodeMap<LimitValueType> b_i(g);
      {
         size_t idx = 0;
         for (int i = 0; i < s; ++i)
            for (int j = 0; j < s; ++j) {
               b_i[U[ID(i, j)]] = +LimitValueType(h1.get(i, j));
               b_i[V[ID(i, j)]] = LimitValueType(0);
               b_i[W[ID(i, j)]] = -LimitValueType(h2.get(i, j));
            }
      }

      // Add all edges
      for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
         const auto& a = arcs[i];
         l_i[a] = 0;
         u_i[a] = simplex.INF;
         c_i[a] = a_costs[i];
      }

      //set lower/upper bounds, cost
      simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

      // Solve the problem to compute the distance
      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType ret = simplex.run();

      switch (ret) {
      case NetworkSimplex<LemonGraph>::INFEASIBLE:
         logger.info("NetworkSimplex<Graph>::INFEASIBLE");
         break;
      case NetworkSimplex<LemonGraph>::OPTIMAL:
         logger.info("NetworkSimplex<Graph>::OPTIMAL");
         break;
      case NetworkSimplex<LemonGraph>::UNBOUNDED:
         logger.info("NetworkSimplex<Graph>::UNBOUNDED");
         break;
      }

      return simplex.totalCost();
   }
   if (h1.getD() == 3) {
      size_t d = s * s * s;

      auto ID = [&s](size_t x, size_t y, size_t z) {
         return x * s * s + y*s + z;
      };

      // add d nodes for each histrogam (d+1) source, (d+2) target
      std::vector<LemonGraph::Node> U;  // left partition
      std::vector<LemonGraph::Node> V;  // middle partition
      std::vector<LemonGraph::Node> W;  // right partition
      std::vector<LemonGraph::Node> Z;  // right partition

      // add first d source nodes
      for (size_t i = 0; i < d; ++i)
         U.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         V.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         W.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         Z.emplace_back(g.addNode());

      std::vector<LemonGraph::Arc> arcs;
      std::vector<int64_t> a_costs;

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int h = 0; h < s; ++h)
               for (int v = 0; v < s; ++v) {
                  arcs.emplace_back(g.addArc(U[ID(i, j, h)], V[ID(i, j, v)]));
                  a_costs.emplace_back(static_cast<int64_t>(pow(h - v, 2)));
               }

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int h = 0; h < s; ++h)
               for (int v = 0; v < s; ++v) {
                  arcs.emplace_back(g.addArc(V[ID(i, j, h)], W[ID(i, v, h)]));
                  a_costs.emplace_back(static_cast<int64_t>(pow(j - v, 2)));
               }

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int h = 0; h < s; ++h)
               for (int v = 0; v < s; ++v) {
                  arcs.emplace_back(g.addArc(W[ID(i, j, h)], Z[ID(v, j, h)]));
                  a_costs.emplace_back(static_cast<int64_t>(pow(i - v, 2)));
               }

      logger.info("Input graph created with %d nodes and %d arcs", countNodes(g), countArcs(g));

      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);

      // lower and upper bounds, cost
      ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
      ListDigraph::ArcMap<LimitValueType> c_i(g);

      // FLow balance
      ListDigraph::NodeMap<LimitValueType> b_i(g);
      {
         size_t idx = 0;
         for (int i = 0; i < s; ++i)
            for (int j = 0; j < s; ++j)
               for (int h = 0; h < s; ++h) {
                  b_i[U[ID(i, j, h)]] = +LimitValueType(h1.get(i, j, h));
                  b_i[V[ID(i, j, h)]] = +LimitValueType(0);
                  b_i[W[ID(i, j, h)]] = -LimitValueType(0);
                  b_i[Z[ID(i, j, h)]] = -LimitValueType(h2.get(i, j, h));
               }
      }

      // Add all edges
      for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
         const auto& a = arcs[i];
         l_i[a] = 0;
         u_i[a] = simplex.INF;
         c_i[a] = a_costs[i];
      }

      //set lower/upper bounds, cost
      simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

      // Solve the problem to compute the distance
      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType ret = simplex.run();

      switch (ret) {
      case NetworkSimplex<LemonGraph>::INFEASIBLE:
         logger.info("NetworkSimplex<Graph>::INFEASIBLE");
         break;
      case NetworkSimplex<LemonGraph>::OPTIMAL:
         logger.info("NetworkSimplex<Graph>::OPTIMAL");
         break;
      case NetworkSimplex<LemonGraph>::UNBOUNDED:
         logger.info("NetworkSimplex<Graph>::UNBOUNDED");
         break;
      }

      return simplex.totalCost();
   }

   if (h1.getD() == 4) {
      size_t d = s * s * s * s;

      auto ID = [&s](size_t x, size_t y, size_t z, size_t u) {
         return x * s * s * s+ y * s * s + z * s + u;
      };

      // add d nodes for each histrogam (d+1) source, (d+2) target
      std::vector<LemonGraph::Node> U;  // left partition
      std::vector<LemonGraph::Node> V;  // middle partition
      std::vector<LemonGraph::Node> W;  // right partition
      std::vector<LemonGraph::Node> Z;  // right partition
      std::vector<LemonGraph::Node> K;  // right partition

      // add first d source nodes
      for (size_t i = 0; i < d; ++i)
         U.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         V.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         W.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         Z.emplace_back(g.addNode());
      for (size_t i = 0; i < d; ++i)
         K.emplace_back(g.addNode());

      std::vector<LemonGraph::Arc> arcs;
      std::vector<int64_t> a_costs;

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int h = 0; h < s; ++h)
               for (int l = 0; l < s; ++l)
                  for (int v = 0; v < s; ++v) {
                     arcs.emplace_back(g.addArc(U[ID(i, j, h, l)], V[ID(i, j, h, v)]));
                     a_costs.emplace_back(static_cast<int64_t>(pow(l - v, 2)));
                  }

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int h = 0; h < s; ++h)
               for (int l = 0; l < s; ++l)
                  for (int v = 0; v < s; ++v) {
                     arcs.emplace_back(g.addArc(V[ID(i, j, h, l)], W[ID(i, j, v, l)]));
                     a_costs.emplace_back(static_cast<int64_t>(pow(h - v, 2)));
                  }

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int h = 0; h < s; ++h)
               for (int l = 0; l < s; ++l)
                  for (int v = 0; v < s; ++v) {
                     arcs.emplace_back(g.addArc(W[ID(i, j, h, l)], Z[ID(i, v, h, l)]));
                     a_costs.emplace_back(static_cast<int64_t>(pow(j - v, 2)));
                  }

      for (int i = 0; i < s; ++i)
         for (int j = 0; j < s; ++j)
            for (int h = 0; h < s; ++h)
               for (int l = 0; l < s; ++l)
                  for (int v = 0; v < s; ++v) {
                     arcs.emplace_back(g.addArc(Z[ID(i, j, h, l)], K[ID(v, j, h, l)]));
                     a_costs.emplace_back(static_cast<int64_t>(pow(i - v, 2)));
                  }

      logger.info("Input graph created with %d nodes and %d arcs", countNodes(g), countArcs(g));

      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> simplex(g);

      // lower and upper bounds, cost
      ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
      ListDigraph::ArcMap<LimitValueType> c_i(g);

      // FLow balance
      ListDigraph::NodeMap<LimitValueType> b_i(g);
      {
         size_t idx = 0;
         for (int i = 0; i < s; ++i)
            for (int j = 0; j < s; ++j)
               for (int h = 0; h < s; ++h)
                  for (int l = 0; l < s; ++l) {
                     b_i[U[ID(i, j, h, l)]] = +LimitValueType(h1.get(i, j, h, l));
                     b_i[V[ID(i, j, h, l)]] = LimitValueType(0);
                     b_i[W[ID(i, j, h, l)]] = LimitValueType(0);
                     b_i[Z[ID(i, j, h, l)]] = LimitValueType(0);
                     b_i[K[ID(i, j, h, l)]] = -LimitValueType(h2.get(i, j, h, l));
                  }
      }

      // Add all edges
      for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
         const auto& a = arcs[i];
         l_i[a] = 0;
         u_i[a] = simplex.INF;
         c_i[a] = a_costs[i];
      }

      //set lower/upper bounds, cost
      simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

      // Solve the problem to compute the distance
      NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>::ProblemType ret = simplex.run();

      switch (ret) {
      case NetworkSimplex<LemonGraph>::INFEASIBLE:
         logger.info("NetworkSimplex<Graph>::INFEASIBLE");
         break;
      case NetworkSimplex<LemonGraph>::OPTIMAL:
         logger.info("NetworkSimplex<Graph>::OPTIMAL");
         break;
      case NetworkSimplex<LemonGraph>::UNBOUNDED:
         logger.info("NetworkSimplex<Graph>::UNBOUNDED");
         break;
      }

      return simplex.totalCost();
   }

}

}