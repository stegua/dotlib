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

namespace DOT {

class SolverNSC {
 public:
  SolverNSC(void) {}

  // Free all pointers
  ~SolverNSC(void) {
    if (!arcs) {
      fprintf(stderr, "KA BOOM!\n");
      return;
    }
    if (!data) {
      fprintf(stderr, "KA BOOM!\n");
      return;
    }
    free((double*)data);
    free((size_t*)arcs);
    arcs = (size_t*)nullptr;
    data = (double*)nullptr;
  }

  // Reserve the internal memory
  void reserveArcs(size_t _n, size_t _m) {
    n = _n;
    m = _m;
    idx = 0;

    size_arcs = 2 * _m;
    arcs = (size_t*)malloc(size_arcs * sizeof(size_t));
    head = &arcs[0];
    tail = &arcs[_m];

    size_data = n + 3 * _m;
    data = (double*)malloc(size_data * sizeof(double));
    cost = &data[0];
    capUB = &data[1 * _m];
    supply = &data[2 * m];

    // Set all suply to zero
    memset(supply, 0, n * sizeof(double));
  }

  // Reset memory for reusing the allocated memory
  void reset(void) {
    memset(arcs, 0, size_arcs * sizeof(size_t));
    memset(data, 0, size_data * sizeof(double));
  }

  // Add node
  void addNode(size_t i, double b) { supply[i] = b; }
  void pluNode(size_t i, double b) { supply[i] += b; }
  void subNode(size_t i, double b) { supply[i] -= b; }

  // Add arc
  void addArc(size_t i, size_t j, double _cost, double _capUB) {
    head[idx] = i;
    tail[idx] = j;

    cost[idx] = _cost;
    capUB[idx] = _capUB;
    idx++;
  }

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
          n = (size_t)std::atol(buf);
          buf = strchr(buf, sep);
          ++buf;
          m = (size_t)std::atol(buf);
          buf = strchr(buf, eol);
          ++buf;
          fprintf(stdout, "%d %d\n", n, m);
          break;
        }
      }
    }

    // check
    if (n > 0 && m > 0) reserveArcs(n, m);

    // Read remaining data
    double tot_remove = 0;
    size_t tot_m = 0;
    offset_cost = 0.0;
    while (buf != NULL) {
      if (buf[0] == 'n') {
        buf = strchr(buf, sep);
        ++buf;

        size_t _idx = (size_t)std::atol(buf);

        buf = strchr(buf, sep);
        ++buf;
        double _s = std::atof(buf);

        buf = strchr(buf, eol);
        if (buf == NULL) break;
        ++buf;

        addNode(_idx - 1, _s);
      } else {
        if (buf[0] == 'a') {
          buf = strchr(buf, sep);
          ++buf;

          size_t i = (size_t)std::atol(buf);

          buf = strchr(buf, sep);
          ++buf;
          size_t j = (size_t)std::atol(buf);

          buf = strchr(buf, sep);
          ++buf;
          double _lb = std::atof(buf);

          buf = strchr(buf, sep);
          ++buf;
          double _ub = std::atof(buf);

          buf = strchr(buf, sep);
          ++buf;
          double _c = std::atof(buf);

          if (_lb == 0.0)
            addArc(i - 1, j - 1, _c, _ub);
          else {
            subNode(i - 1, _lb);
            pluNode(j - 1, _lb);
            addArc(i - 1, j - 1, _c, _ub - _lb);

            offset_cost += _c * _lb;
          }

          tot_m++;
          if (tot_m == m) break;

          buf = strchr(buf, eol);
          ++buf;
        }
      }
    }

    return;
  }

  // Solve the instance
  double solve() {
    m = idx;
    NetSimplexCapacity simplex('F', n, m);
    simplex.setTimelimit(3600);
    simplex.setOptTolerance(0.0);

    for (size_t i = 0; i < n; i++) simplex.addNode((int)i, supply[i]);

    for (size_t e = 0; e < m; e++)
      simplex.addArc((int)head[e], (int)tail[e], cost[e], capUB[e], e);

    simplex.run();

    double r = offset_cost + simplex.totalCost();

    size_t api_it = simplex.iterations();
    double api_time = simplex.runtime();

    fprintf(stdout, "CoinLEM %.f ==> %d#%d, it: %d, time: %.3f\n", r, n, m,
            api_it, api_time);

    return r;
  }

  // solve the instance by column generation
  double colgen() {
    auto start_t = std::chrono::steady_clock::now();
    m = idx;

    std::unordered_set<int> Ks;
    for (int e = 0; e < m; ++e) Ks.insert(tail[e]);
    int nk = Ks.size();
    unordered_map<int, int> M;
    int vv = 0;
    for (int u : Ks) {
      M[u] = vv;
      ++vv;
    }

    CapVars vnew;
    vnew.reserve(nk);

    vector<double> best_v(nk, 0);
    vector<int> best_e(nk, -1);

    // Dual multipliers
    vector<double> pi(n, 0);

    // Simplex
    NetSimplexCapacity simplex('E', n, 0);
    simplex.setTimelimit(3600);

    double negeps = std::nextafter(-1e-09, -0.0);
    simplex.setOptTolerance(negeps);

    double ART_COST = 0;
    for (int e = 0; e < m; ++e)
      if (cost[e] > ART_COST) ART_COST = cost[e];
    ART_COST = (ART_COST + 1) * n;
    simplex.setArtCost(ART_COST);

    for (size_t i = 0; i < n; i++) simplex.addNode((int)i, supply[i]);

    ProblemType _status = simplex.run();

    std::vector<bool> ARCS(m, false);

    double _all_p = 0.0;
    while (_status != ProblemType::TIMELIMIT) {
      // Take the dual values
      for (int j = 0; j < n; ++j) pi[j] = simplex.potential(j);

      for (int k = 0; k < nk; ++k) {
        best_v[k] = negeps;
        best_e[k] = -1;
      }

      // Solve separation problem:
      auto start_p = std::chrono::steady_clock::now();

      // To improve this loop ----------------------
      for (int e = 0; e < m; ++e) {
        int k = M[tail[e]];
        double violation = cost[e] + pi[head[e]] - pi[tail[e]];
        if (violation < best_v[k] && !ARCS[e]) {
          best_v[k] = violation;
          best_e[k] = e;
        }
      }

      auto end_p = std::chrono::steady_clock::now();
      _all_p += double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_p - start_p)
                           .count()) /
                1000;
      //---------------------------------------------

      // Take all negative reduced cost variables
      vnew.clear();
      for (int k = 0; k < nk; ++k)
        if (best_e[k] >= 0) {
          int e = best_e[k];
          vnew.emplace_back(head[e], tail[e], cost[e], capUB[e], e);
          ARCS[e] = true;
        }

      if (vnew.empty()) {
        fprintf(stdout, " no new vars\n");
        break;
      }

      // std::sort(vnew.begin(), vnew.end(),
      //          [](const CapVar& v, const CapVar& w) { return v.c > w.c; });

      // Replace old constraints with new ones
      int new_arcs = simplex.addArcs(vnew);
      if (new_arcs == 0) break;

      _status = simplex.reRun();
    }

    double r = offset_cost + simplex.totalCost();

    size_t api_it = simplex.iterations();
    double api_time = simplex.runtime();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    fprintf(stdout,
            "CHECK SOL %.3f ==> %d#%d, it: %d, alltime: %.3f, simplex: %.3f, "
            "all_p: %.3f, arcs: %d\n",
            r, n, m, api_it, _all, api_time, _all_p, simplex.num_arcs());

    return r;
  }

  // solve the instance by column generation
  double block_colgen() {
    auto start_t = std::chrono::steady_clock::now();

    m = idx;

    int NUM_BLOCKS = 1024;
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
    simplex.setOptTolerance(negeps);

    double ART_COST = 0;
    for (int e = 0; e < m; ++e)
      if (cost[e] > ART_COST) ART_COST = cost[e];
    ART_COST = (ART_COST + 1) * n;
    simplex.setArtCost(ART_COST);

    for (size_t i = 0; i < n; i++) simplex.addNode((int)i, supply[i]);

    std::vector<bool> ARCS(m, false);

    //// Heuristic columns
    // vector<int> start_pivots;
    // start_pivots.reserve(n);
    // for (size_t i = 0; i < n; i++)
    //  if (supply[i] > 0) {
    //    double c, min_cost = std::numeric_limits<double>::max();
    //    int min_arc = -1;
    //    for (int e = 0; e < m; ++e)
    //      if (head[e] == i) {
    //        c = cost[e];
    //        if (c < min_cost) {
    //          min_cost = c;
    //          min_arc = e;
    //        }
    //      }
    //    if (min_arc > -1) start_pivots.push_back(min_arc);
    //  }

    //// Initial heuristics
    // fprintf(stdout, "start %d\n", (int)start_pivots.size());
    // for (int e : start_pivots) {
    //  simplex.addArc(head[e], tail[e], cost[e], capUB[e], e);
    //  ARCS[e] = true;
    //}

    ProblemType _status = simplex.run();

    double _all_p = 0.0;

    while (_status != ProblemType::TIMELIMIT) {
      // Solve separation problem:

      // Take the dual values
      for (int j = 0; j < n; ++j) pi[j] = simplex.potential(j);

      // To improve this loop ----------------------
      // auto start_p = std::chrono::steady_clock::now();
      double itime = omp_get_wtime();
      //#pragma omp parallel for
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
              break;
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
          vnew.emplace_back(head[e], tail[e], cost[e], capUB[e], e);
          ARCS[e] = true;
        }

      if (vnew.empty()) {
        fprintf(stdout, " no new vars. block\n");
        break;
      }

      // std::sort(vnew.begin(), vnew.end(),
      //          [](const CapVar& v, const CapVar& w) { return v.c < w.c; });

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
            "CHECK SOL %.3f ==> %d#%d, it: %d, alltime: %.3f, simplex: %.3f, "
            "all_p: %.3f, arcs: %d\n",
            r, n, m, api_it, _all, api_time, _all_p, simplex.num_arcs());

    return r;
  }

  // Dump the graph to a file
  void dump(void) const {
    for (size_t u = 0; u < n; u++) fprintf(stdout, "%ld %f\n", u, supply[u]);

    for (size_t e = 0; e < m; e++)
      fprintf(stdout,
              "(head: %ld, tail: %ld, lb: %.3f, ub: %.3f, cost: %.3f)\n",
              head[e], tail[e], capUB[e], cost[e]);
  }

  // Class attributes
  size_t n;
  size_t m;

  size_t idx;

  // Keep contigous memory
  size_t size_arcs;
  size_t* arcs;
  size_t* head;
  size_t* tail;

  // Offset of cost for handling lower bounds
  double offset_cost;

  // Keep contigous memory
  size_t size_data;
  double* data;
  double* cost;    // one for arcs
  double* capUB;   // one for arcs
  double* supply;  // one for node
};                 // namespace DOT
}  // namespace DOT
