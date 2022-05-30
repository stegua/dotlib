/**
 * @fileoverview Copyright (c) 2019-2022, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

// ORIGINAL SOURCE CODE FOR THE NETWORK SIMPLEX BASIS DATA STRUCTURE TAKE FROM:
// WEBSITE: https://lemon.cs.elte.hu

/* ORIGINAL LICENSE FILE:
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2013
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#pragma once

#include "DOT_Commons.h"

// Cplex Network Simplex
#define LEMON_ONLY_TEMPLATES
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>

typedef lemon::ListDigraph LemonGraph;
typedef double LimitValueType;
typedef lemon::NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>
    LemonSimplex;

namespace DOT {

class SolverCoinLemon {
 public:
  SolverCoinLemon(void) {}

  // Free all pointers
  ~SolverCoinLemon(void) {}

  // Parse a DIMACS file
  void parseDIMACS(const std::string filename, size_t& len) {
    /* EXAMPLE of DIMACS file from Austin:
       p min 86 1223
       n 1 13002700
       n 86 -13002700
       a 1 2 30990 30990 0
       a 1 3 64749 64749 0
       a 1 4 122822 122822 0
       a 1 5 74129 74129 0
    */

    const char sep = ' ';
    const char eol = '\n';

    std::string raw = readTextFile(filename, len);

    const char* buf = raw.data();

    // read graph dimensions
    while (buf != NULL) {
      if (buf[0] == 'c') {
        buf = strchr(buf, sep);
        buf = strchr(buf, eol);
      } else {
        if (buf[0] == 'p') {
          buf = strchr(buf, sep);
          ++buf;
          buf = strchr(buf, sep);
          ++buf;
          n = std::atoll(buf);
          buf = strchr(buf, sep);
          ++buf;
          m = std::atoll(buf);
          buf = strchr(buf, eol);
          ++buf;
          fprintf(stdout, "%d %d\n", n, m);
          break;
        }
      }
    }

    // check
    // if (n > 0 && m > 0) reserveArcs(n, m);
    G.reserveNode(n);
    supply.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      nodes.emplace_back(G.addNode());
      supply.emplace_back(0.0);
    }

    G.reserveArc(m);
    arcs.reserve(m);
    cost.reserve(m);
    lb.reserve(m);
    ub.reserve(m);

    // Read remaining data

    while (buf[0] == 'n') {
      buf = buf + 2;

      size_t _idx = std::atoll(buf);

      buf = std::strchr(buf, sep);
      ++buf;
      double _s = std::atof(buf);

      buf = std::strchr(buf, eol);
      if (buf == NULL) break;
      ++buf;

      supply[_idx - 1] = _s;
    }

    size_t i, j;
    double _c, _lb, _ub;

    while (buf[0] == 'a') {
      buf = buf + 2;

      i = std::atoll(buf);

      buf = std::strchr(buf, sep);
      ++buf;
      j = std::atoll(buf);

      buf = std::strchr(buf, sep);
      ++buf;
      _lb = std::atof(buf);

      buf = std::strchr(buf, sep);
      ++buf;
      _ub = std::atof(buf);

      buf = std::strchr(buf, sep);
      ++buf;
      _c = std::atof(buf);
      arcs.emplace_back(G.addArc(nodes[i - 1], nodes[j - 1]));
      cost.emplace_back(_c);
      lb.emplace_back(_lb);
      ub.emplace_back(_ub);

      buf = std::strchr(buf, eol);
      ++buf;
    }
  }

  // Solve the instance
  double solve(const std::string& msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    LemonSimplex simplex(G, false);

    // lower and upper bounds, cost
    LemonGraph::ArcMap<LimitValueType> l_i(G), u_i(G);
    LemonGraph::ArcMap<LimitValueType> c_i(G);

    // FLow balance
    LemonGraph::NodeMap<LimitValueType> b_i(G);

    for (size_t i = 0; i < n; ++i) b_i[nodes[i]] = supply[i];

    // Add all edges
    for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = lb[i];
      u_i[a] = ub[i];
      c_i[a] = cost[i];
    }

    // set lower/upper bounds, cost
    simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

    // Solve the problem to compute the distance
    lemon::NetworkSimplex<LemonGraph, LimitValueType,
                          LimitValueType>::ProblemType ret = simplex.run();

    switch (ret) {
      case lemon::NetworkSimplex<LemonGraph>::INFEASIBLE:
        fprintf(stdout, "NetworkSimplex<Graph>::INFEASIBLE\n");
        break;
      case lemon::NetworkSimplex<LemonGraph>::OPTIMAL:
        fprintf(stdout, "NetworkSimplex<Graph>::OPTIMAL\n");
        break;
      case lemon::NetworkSimplex<LemonGraph>::UNBOUNDED:
        fprintf(stdout, "NetworkSimplex<Graph>::UNBOUNDED\n");
        break;
    }

    double objval = simplex.totalCost();

    // Merge with previous code
    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    double _iterations = 0;

    fprintf(stdout,
            "LEM-PLEX %s it %d UB %.6f runtime %.4f simplex %.4f "
            "num_arcs %d\n",
            msg.c_str(), (int)_iterations, objval, _all, _all, (int)m);

    return _all;
  }

  // Dump the graph to a file
  void dump(void) const {
    for (size_t u = 0; u < n; u++) fprintf(stdout, "%ld %f\n", u, supply[u]);

    // for (size_t e = 0; e < m; e++)
    //  fprintf(stdout,
    //          "(head: %ld, tail: %ld, lb: %.3f, ub: %.3f, cost: %.3f)\n",
    //          (int)head[e], (int)tail[e], capLB[e], capUB[e], cost[e]);
  }

  // Class attributes
  size_t n;
  size_t m;

  LemonGraph G;
  std::vector<LemonGraph::Node> nodes;
  std::vector<LemonGraph::Arc> arcs;

  vector<LimitValueType> supply;
  vector<double> cost;
  vector<double> lb;
  vector<double> ub;
};
}  // namespace DOT
