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

#include <omp.h>

#include "DOT_Commons.h"
#include "DOT_NetSimplexCap.h"

// Cplex Network Simplex
#define LEMON_ONLY_TEMPLATES
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>

typedef lemon::ListDigraph LemonGraph;
typedef double LimitValueType;
typedef lemon::NetworkSimplex<LemonGraph, LimitValueType, LimitValueType>
    LemonSimplex;

namespace DOT {

class SolverTRI {
 public:
  SolverTRI(void) {}

  // Free all pointers
  ~SolverTRI(void) {}

  // Parse a DIMACS file
  void parseDIMACS(const std::string filename, size_t &len) {
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

    const char *buf = raw.data();

    offset_cost = 0.0;

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
      ub.emplace_back(_ub - _lb);

      supply[i - 1] -= _lb;
      supply[j - 1] += _lb;

      offset_cost += _c * _lb;

      buf = std::strchr(buf, eol);
      ++buf;
    }
  }

  // Parse a DIMACS file
  void parseDIMACSTRI(const std::string filename, size_t &len) {
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

    const char *buf = raw.data();

    offset_cost = 0.0;

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
          n = std::atoi(buf);
          buf = strchr(buf, sep);
          ++buf;
          m = std::atoi(buf);
          buf = strchr(buf, sep);
          ++buf;
          k = std::atoi(buf);
          buf = strchr(buf, eol);
          ++buf;
          fprintf(stdout, "%d %d %d\n", n, m, k);
          break;
        }
      }
    }

    // check
    // if (n > 0 && m > 0) reserveArcs(n, m);
    supply.resize(n, 0.0);

    cost.resize(m, 0.0);
    ub.resize(m, 0.0);

    head.resize(m, -1);
    tail.resize(m, -1);

    // Read remaining data

    while (buf[0] == 'n') {
      buf = buf + 2;

      int _idx = std::atoi(buf);

      buf = std::strchr(buf, sep);
      ++buf;
      double _s = std::atof(buf);

      buf = std::strchr(buf, eol);
      if (buf == NULL) break;
      ++buf;

      if (_idx == n) _idx = 0;
      supply[_idx] = _s;
    }

    size_t i, j;
    double _c, _lb, _ub;

    int e = 0;
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

      if (j == n) {
        e = i - j + k;
        j = 0;

        head[e] = i;
        tail[e] = j;
        cost[e] = _c;
        ub[e] = _ub - _lb;

        supply[i] -= _lb;
        supply[j] += _lb;
      } else {
        e = i * k + (j - (n - 1 - k)) - 1;

        head[e] = i;
        tail[e] = j;
        cost[e] = _c;
        ub[e] = _ub - _lb;

        supply[i] -= _lb;
        supply[j] += _lb;
      }

      offset_cost += _c * _lb;

      buf = std::strchr(buf, eol);
      ++buf;
    }

    // Sort arcs by cost per blocks
    vector<int> Es;
    Es.reserve(m);
    for (int e = 0; e < m; ++e) Es.push_back(e);

    std::random_device rd;
    std::mt19937 g(13);

    // std::shuffle(Es.begin() + k, Es.end(), g);

    // for (int e = k; e < m; e = e + k) {
    //  std::sort(&Es[e], &Es[e + k],
    //            [&](const int &v, const int &w) { return cost[v] > cost[w];
    //            });
    //}
    //  for (int e = 0; e < m; ++e) {
    //  fprintf(stdout, "%d %.f\n", e, cost[e]);
    //}
    // fprintf(stdout, "\n");
    // for (int e = 0; e < m; ++e) {
    //  fprintf(stdout, "%d %.f\n", Es[e], cost[Es[e]]);
    //}
    vector<double> tcost(cost);
    vector<double> tub(ub);
    vector<int> thead(head);
    vector<int> ttail(tail);
    for (int e = k; e < m; ++e) {
      tcost[e] = cost[Es[e]];
      tub[e] = ub[Es[e]];
      thead[e] = head[Es[e]];
      ttail[e] = tail[Es[e]];
    }

    std::swap(ub, tub);
    std::swap(cost, tcost);
    std::swap(head, thead);
    std::swap(tail, ttail);
  }

  // Solve the instance
  double solve(const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    NetSimplexCapacity simplex('F', n, m);
    simplex.setTimelimit(3600);

    simplex.setOptTolerance(1e-09);

    double ART_COST = 0;
    for (int e = 0; e < m; ++e)
      if (cost[e] > ART_COST) ART_COST = cost[e];
    ART_COST = (ART_COST + 1) * n;
    simplex.setArtCost(ART_COST);

    for (int i = 0; i < n; i++) {
      simplex.addNode((int)i, supply[i]);
      // fprintf(stdout, "%.f\n", supply[i]);
    }
    // Store the arcs in the original order
    for (int e = 0; e < m; ++e) {
      simplex.addArc(head[e], tail[e], cost[e], ub[e], e);
      // fprintf(stdout, "id: %d -> %d %d %.f %.f\n", e, head[e], tail[e],
      // cost[e],
      //        ub[e]);
    }

    simplex.run();

    double objval = offset_cost + simplex.totalCost();

    size_t api_it = simplex.iterations();
    double api_time = simplex.runtime();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    fprintf(stdout,
            "TRI-PLEX %s it %d UB %.f runtime %.4f simplex %.4f "
            "num_arcs %d\n",
            msg.c_str(), (int)api_it, objval, _all, api_time, (int)m);

    return _all;
  }

  // solve the instance by column generation
  double colgen() {
    auto start_t = std::chrono::steady_clock::now();

    int NUM_BLOCKS = n - k - 1;
    int block = m / NUM_BLOCKS + (m % NUM_BLOCKS != 0);
    fprintf(stdout, "m: %d, blocks: %d, size_blocks: %d\n", m, NUM_BLOCKS,
            block);

    CapVars vnew;
    vnew.reserve(NUM_BLOCKS);

    vector<double> best_v(NUM_BLOCKS, 0);
    vector<int> best_e(NUM_BLOCKS, -1);

    // Dual multipliers
    vector<double> pi(n, 0);

    // Simplex
    NetSimplexCapacity simplex('E', n, 0);
    simplex.setTimelimit(3600);

    double negeps = std::nextafter(-1e-09, -0.0);
    simplex.setOptTolerance(1e-09);

    double ART_COST = 0;
    for (int e = 0; e < m; ++e)
      if (cost[e] > ART_COST) ART_COST = cost[e];
    ART_COST = (ART_COST + 1) * n;
    simplex.setArtCost(ART_COST);

    for (size_t i = 0; i < n; i++) simplex.addNode((int)i, supply[i]);

    std::vector<bool> ARCS(m, false);

    for (int e = 0; e < k; ++e) {
      simplex.addArc(head[e], tail[e], cost[e], ub[e], e);
      ARCS[e] = true;
    }

    ProblemType _status = simplex.run();

    double _all_p = 0.0;

    while (_status != ProblemType::TIMELIMIT) {
      // Solve separation problem:

      // Take the dual values
      for (int j = 0; j < n; ++j) pi[j] = simplex.potential(j);

      // To improve this loop ----------------------
      // auto start_p = std::chrono::steady_clock::now();
      double itime = omp_get_wtime();
      //#pragma omp parallel for schedule(static, 1) num_threads(8)
      for (int b = 0; b < NUM_BLOCKS; ++b) {
        best_v[b] = negeps;
        best_e[b] = -1;

        for (int e = b * block, e_max = std::min<int>(m, (b + 1) * block);
             e < e_max; ++e) {
          if (!ARCS[e]) {
            double violation = cost[e] + pi[head[e]] - pi[tail[e]];
            if (violation < best_v[b]) {
              best_v[b] = violation;
              best_e[b] = e;
            }
          }
        }
      }
      // auto end_p = std::chrono::steady_clock::now();
      _all_p += (omp_get_wtime() - itime);
      /*double(std::chrono::duration_cast<std::chrono::milliseconds>(
                          end_p - start_p)
                          .count()) /
               1000;*/

      //---------------------------------------------

      // Take all negative reduced cost variables
      vnew.clear();
      for (int b = 0; b < NUM_BLOCKS; ++b)
        if (best_e[b] > -1) {
          int e = best_e[b];
          vnew.emplace_back(head[e], tail[e], cost[e], ub[e], e);
          ARCS[e] = true;
        }

      if (vnew.empty()) {
        fprintf(stdout, " no new vars. block\n");
        break;
      }

      // std::sort(vnew.begin(), vnew.end(),
      //          [](const CapVar &v, const CapVar &w) { return v.c < w.c; });

      // Replace old constraints with new ones
      simplex.updateArcs(vnew, ARCS);

      _status = simplex.reRun();
    }

    double r = simplex.totalCost();

    size_t api_it = simplex.iterations();
    double api_time = simplex.runtime();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    fprintf(stdout,
            "BLOCK-CG %.f ==> %d#%d, it: %d, alltime: %.3f, simplex: %.3f, "
            "all_p: %.3f, arcs: %d\n",
            r, n, m, api_it, _all, api_time, _all_p, simplex.num_arcs());

    return r;
  }

  double block_colgen() {
    auto start_t = std::chrono::steady_clock::now();

    int NUM_BLOCKS = n - k;
    int block = k;
    fprintf(stdout, "m: %d, blocks: %d, size_blocks: %d\n", m, NUM_BLOCKS,
            block);

    CapVars vnew;
    vnew.reserve(NUM_BLOCKS);

    vector<double> best_v(NUM_BLOCKS, 0);
    vector<int> best_e(NUM_BLOCKS, -1);

    // Dual multipliers
    vector<double> pi(n, 0);

    // Simplex
    NetSimplexCapacity simplex('E', n, 0);
    simplex.setTimelimit(3600);

    double negeps = std::nextafter(-1e-09, -0.0);
    simplex.setOptTolerance(1e-09);

    double ART_COST = 0;
    for (int e = 0; e < m; ++e)
      if (cost[e] > ART_COST) ART_COST = cost[e];
    ART_COST = (ART_COST + 1) * n;
    simplex.setArtCost(ART_COST);

    for (size_t i = 0; i < n; i++) simplex.addNode((int)i, supply[i]);

    std::vector<bool> ARCS(m, false);

    for (int e = 0; e < k; ++e) {
      simplex.addArc(head[e], tail[e], cost[e], ub[e], e);
      ARCS[e] = true;
    }

    ProblemType _status = simplex.run();

    double _all_p = 0.0;

    while (_status != ProblemType::TIMELIMIT) {
      // Solve separation problem:

      // Take the dual values
      for (int j = 0; j < n; ++j) pi[j] = simplex.potential(j);

      // To improve this loop ----------------------
      // auto start_p = std::chrono::steady_clock::now();
      double itime = omp_get_wtime();

      for (int b = 0; b < NUM_BLOCKS; ++b) {
        best_v[b] = negeps;
        best_e[b] = -1;

        for (int l = 0; l < k; ++l) {
          int e = b * k + l;
          if (!ARCS[e]) {
            double violation = cost[e] + pi[head[e]] - pi[tail[e]];
            if (violation < best_v[b]) {
              best_v[b] = violation;
              best_e[b] = e;
            }
          }
        }
      }
      // auto end_p = std::chrono::steady_clock::now();
      _all_p += (omp_get_wtime() - itime);
      /*double(std::chrono::duration_cast<std::chrono::milliseconds>(
                          end_p - start_p)
                          .count()) /
               1000;*/

      //---------------------------------------------

      // Take all negative reduced cost variables
      vnew.clear();
      for (int b = 0; b < NUM_BLOCKS; ++b)
        if (best_e[b] > -1) {
          int e = best_e[b];
          vnew.emplace_back(head[e], tail[e], cost[e], ub[e], e);
          ARCS[e] = true;
        }

      if (vnew.empty()) {
        fprintf(stdout, " no new vars. block\n");
        break;
      }

      // std::sort(vnew.begin(), vnew.end(),
      //          [](const CapVar &v, const CapVar &w) { return v.c < w.c; });

      // Replace old constraints with new ones
      simplex.updateArcs(vnew, ARCS);

      _status = simplex.reRun();
    }

    double r = simplex.totalCost();

    size_t api_it = simplex.iterations();
    double api_time = simplex.runtime();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    fprintf(stdout,
            "BLOCK-CG %.f ==> %d#%d, it: %d, alltime: %.3f, simplex: %.3f, "
            "all_p: %.3f, arcs: %d\n",
            r, n, m, api_it, _all, api_time, _all_p, simplex.num_arcs());

    return r;
  }

  // Class attributes
  int n;
  int m;
  int k;

  // Offset of cost for handling lower bounds
  double offset_cost;

  // Input graph
  LemonGraph G;
  std::vector<LemonGraph::Node> nodes;
  std::vector<LemonGraph::Arc> arcs;

  vector<LimitValueType> supply;
  vector<double> cost;
  vector<double> ub;
  vector<int> head;
  vector<int> tail;

  // Trigraph
  vector<int> t_supply;

  vector<int> t_arcs;  // flat matrix
};

}  // namespace DOT
