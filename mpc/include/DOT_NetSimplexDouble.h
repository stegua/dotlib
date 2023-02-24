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

namespace DOT {
static void eati_search_dbl(int e, int e_max, const int *state,
                            const double *cost, const double *pi,
                            const int *source, const int *target, double *min,
                            int *in_arc) {
  // min[0] = 0;
  for (int i = e; i < e_max; ++i) {
    int c = double(state[i]) * (cost[i] + pi[source[i]] - pi[target[i]]);
    if (c < min[0]) {
      min[0] = c;
      in_arc[0] = i;
    }
  }
}

class NetSimplexDouble {
 public:
  typedef std::vector<int> IntVector;
  typedef std::vector<double> ValueVector;
  typedef std::vector<double> CostVector;
  typedef std::vector<signed char> CharVector;

  typedef std::vector<bool> BoolVector;

  // State constants for arcs
  const bool STATE_TREE = false;
  const bool STATE_LOWER = true;

  // Direction constants for tree arcs
  const int DIR_DOWN = -1;
  const int DIR_UP = 1;

  // Data related to the underlying digraph
  int _node_num;
  int _arc_num;

  int _dummy_arc;  // Arc id where begin the basic arcs
  int _next_arc;

  // Parameters of the problem
  double _sum_supply;

  // Data structures for storing the digraph
  IntVector _source;
  IntVector _target;

  // Node and arc data
  ValueVector _supply;
  ValueVector _flow;
  CostVector _cost;
  CostVector _pi;

  // Data for storing the spanning tree structure
  IntVector _parent;
  IntVector _pred;
  IntVector _thread;
  IntVector _rev_thread;
  IntVector _succ_num;
  IntVector _last_succ;
  CharVector _pred_dir;
  IntVector _state;
  IntVector _dirty_revs;

  int _root;

  // Temporary data used in the current pivot iteration
  int in_arc, join, u_in, v_in, u_out, v_out;
  double delta;

  const double MAX;
  const double INF;

  double _runtime;
  double _time_pricing;
  double _time_update_duals;
  double _time_update_basis;

  double _timelimit;
  std::string _verbosity;
  double _opt_tolerance;

  int N_IT_LOG;

  uint64_t _iterations;

  double t1, t2, t3, t4, t5, t6;

  const int BLOCK_SIZE_FACTOR = 4;
  const int MIN_BLOCK_SIZE = 20;
  int _block_size;

 private:
  // Implementation of the Block Search pivot rule
  class BlockSearchPivotRule {
   private:
    // References to the NetworkSimplex class
    const IntVector &_source;
    const IntVector &_target;
    const CostVector &_cost;
    const IntVector &_state;
    const CostVector &_pi;
    int &_in_arc;
    int _arc_num;
    int _dummy_arc;

    // Pivot rule data
    int _block_size;
    int _next_arc;

    // Negative eps
    const double negeps;

   public:
    // Constructor
    BlockSearchPivotRule(NetSimplexDouble &ns)
        : _source(ns._source),
          _target(ns._target),
          _cost(ns._cost),
          _state(ns._state),
          _pi(ns._pi),
          _in_arc(ns.in_arc),
          _arc_num(ns._arc_num),
          _dummy_arc(ns._dummy_arc),
          _next_arc(ns._next_arc),
          negeps(-ns._opt_tolerance) {
      // The main parameters of the pivot rule
      const double BLOCK_SIZE_FACTOR = 4;
      const int MIN_BLOCK_SIZE = 20;

      _block_size =
          (std::max)(int(BLOCK_SIZE_FACTOR *
                         std::sqrt(double(_arc_num) - double(_dummy_arc))),
                     int(BLOCK_SIZE_FACTOR *
                         std::sqrt(double(_source.size()))));
    }

    // Find next entering arc
    bool findEnteringArc() {
      double min = negeps;
      int e = _next_arc;
      int M = _arc_num;
      int M0 = _dummy_arc;
      int e_max = std::min(M, e + _block_size);

      while (e < M) {
        eati_search_dbl(e, e_max, &_state[0], &_cost[0], &_pi[0], &_source[0],
                        &_target[0], &min, &_in_arc);

        if (min < negeps) {
          _next_arc = _in_arc;
          return true;
        }

        e = e_max;
        e_max = std::min(M, e + _block_size);
      }

      e = M0;
      e_max = std::min(_next_arc, e + _block_size);

      while (e < _next_arc) {
        eati_search_dbl(e, e_max, &_state[0], &_cost[0], &_pi[0], &_source[0],
                        &_target[0], &min, &_in_arc);

        if (min < negeps) {
          _next_arc = _in_arc;
          return true;
        }

        e = e_max;
        e_max = std::min(_next_arc, e + _block_size);
      }

      return false;
    }

  };  // class BlockSearchPivotRule

 public:
  NetSimplexDouble(const char INIT, int node_num, int arc_num)
      : _node_num(node_num),
        _arc_num(0),
        _root(-1),
        in_arc(-1),
        join(-1),
        u_in(-1),
        v_in(-1),
        u_out(-1),
        v_out(-1),
        MAX((std::numeric_limits<double>::max)()),
        INF(std::numeric_limits<double>::has_infinity
                ? std::numeric_limits<double>::infinity()
                : MAX),
        _runtime(0.0) {
    // Check the number types
    if (!std::numeric_limits<int>::is_signed)
      throw std::runtime_error(
          "The flow type of NetworkSimplex must be signed");
    if (!std::numeric_limits<int>::is_signed)
      throw std::runtime_error(
          "The cost type of NetworkSimplex must be signed");

    // Reset data structures
    int all_node_num = _node_num + 1;
    _supply.resize(all_node_num, 0);
    _pi.resize(all_node_num);

    _parent.resize(all_node_num);    // Parent
    _pred.resize(all_node_num);      // Predecessor N -> E
    _pred_dir.resize(all_node_num);  // Predecessor direction N -> {-1,1}
    _thread.resize(all_node_num);    // Depth-first order for visiting the tree
    _rev_thread.resize(all_node_num);  // Reverse thread (?)
    _succ_num.resize(all_node_num);    // Number of successors?
    _last_succ.resize(all_node_num);   // Last successor?

    // 2*n arcs from nodes to root and from root to node;
    // 2*n-1 nodes in a basic solution
    int max_arc_num = 0;
    if (INIT == 'F')  // Full
      max_arc_num = _node_num + arc_num;

    if (INIT == 'E')  // Empty, for Column Generation
      max_arc_num = 4 * _node_num + 1;

    _source.reserve(max_arc_num);
    _target.reserve(max_arc_num);

    _cost.reserve(max_arc_num);
    _flow.reserve(max_arc_num);
    _state.reserve(max_arc_num);

    _source.resize(_node_num);
    _target.resize(_node_num);

    _cost.resize(_node_num, 0.0);
    _flow.resize(_node_num, 0.0);
    _state.resize(_node_num, STATE_LOWER);

    // Dummy arcs for eavery node to root node
    _dummy_arc = _node_num;
    _arc_num = _dummy_arc;
    _next_arc = _dummy_arc;

    // Interal parameters
    _timelimit = std::numeric_limits<double>::max();
    _verbosity = DOT_VAL_INFO;
    _opt_tolerance = 1e-09;
    _iterations = 0;

    // Block size for the pivot
    _block_size =
        std::max(int(BLOCK_SIZE_FACTOR * std::sqrt(double(max_arc_num))),
                 MIN_BLOCK_SIZE);
  }

  NetSimplexDouble(const NetSimplexDouble &o)
      : _source(o._source),
        _target(o._target),
        _supply(o._supply.begin(), o._supply.end()),
        _flow(o._flow.begin(), o._flow.end()),
        _cost(o._cost.begin(), o._cost.end()),
        _pi(o._pi.begin(), o._pi.end()),
        _parent(o._parent),
        _pred(o._pred),
        _thread(o._thread),
        _rev_thread(o._rev_thread),
        _succ_num(o._succ_num),
        _last_succ(o._last_succ),
        _pred_dir(o._pred_dir),
        _state(o._state),
        _dirty_revs(o._dirty_revs),
        _dummy_arc(o._dummy_arc),
        _next_arc(o._dummy_arc),
        _node_num(o._node_num),
        _arc_num(o._arc_num),
        _root(o._root),
        in_arc(o.in_arc),
        join(o.join),
        u_in(o.u_in),
        v_in(o.v_in),
        u_out(o.u_out),
        v_out(o.v_out),
        MAX(o.MAX),
        INF(o.INF),
        _runtime(o._runtime) {}

  ProblemType run(PivotRule pivot_rule = PivotRule::BLOCK_SEARCH) {
    // shuffle();
    _runtime = 0.0;
    _iterations = 0;

    // Reset arc variables
    for (int e = 0; e < _arc_num; ++e) {
      _state[e] = STATE_LOWER;
      _flow[e] = 0.0;
    }

    if (!init()) return ProblemType::INFEASIBLE;
    return start(pivot_rule);
  }

  ProblemType reRun(PivotRule pivot_rule = PivotRule::BLOCK_SEARCH) {
    return start(pivot_rule);
  }

  uint64_t num_arcs() const { return uint64_t(_source.size()) - _dummy_arc; }

  uint64_t num_nodes() const { return _node_num; }

  void addNode(int i, double b) { _supply[i] = b; }

  size_t addArc(int a, int b, double c) {
    size_t idx = _source.size();
    _source.emplace_back(a);
    _target.emplace_back(b);
    _cost.emplace_back(c);

    _flow.emplace_back(0);
    _state.emplace_back(STATE_LOWER);

    _arc_num++;
    return idx;
  }

  void setArc(size_t idx, int a, int b, double c) {
    _source[_dummy_arc + idx] = a;
    _target[_dummy_arc + idx] = b;
    _cost[_dummy_arc + idx] = c;

    _flow[_dummy_arc + idx] = 0;
    _state[_dummy_arc + idx] = STATE_LOWER;

    _arc_num++;
  }

  // Manage dummy flow
  void updateArcs(const vector<size_t> &as, double value) {
    for (const auto e : as) _cost[e] = value;
  }

  double computeDummyFlow(const vector<size_t> &as) const {
    double tot_flow = 0;
    for (const auto e : as) tot_flow += _flow[e];
    return tot_flow;
  }

  // Column Generation
  int updateArcs(const DoubleVars &as) {
    int new_arc = 0;
    size_t idx = 0;
    size_t idx_max = as.size();

    int e = _dummy_arc;
    int e_max = _arc_num;

    // Store the new arc variable, replacing an used arc,
    // out of the basis and with positive reduced cost
    for (; idx < idx_max; ++idx) {
      while (e < e_max) {
        // Replace useless variables with new variables
        if (_state[e] == STATE_LOWER &&
            (_cost[e] + _pi[_source[e]] - _pi[_target[e]] > 0))
          break;
        ++e;
      }
      if (e >= e_max) break;
      _source[e] = as[idx].a;
      _target[e] = as[idx].b;
      _cost[e] = as[idx].c;
      if (new_arc == 0) _next_arc = e;
      new_arc++;
    }

    for (; idx < idx_max; ++idx) {
      addArc(as[idx].a, as[idx].b, as[idx].c);
      if (new_arc == 0) _next_arc = e;
      new_arc++;
    }

    return new_arc;
  }

  int addArcs(const DoubleVars &as) {
    int new_arc = 0;
    size_t idx = 0;
    size_t idx_max = as.size();

    for (; idx < idx_max; ++idx) {
      addArc(as[idx].a, as[idx].b, as[idx].c);
      if (new_arc == 0) _next_arc = _arc_num;
      new_arc++;
    }

    return new_arc;
  }

  double totalCost() const {
    double c = 0;
    for (int e = _dummy_arc; e < _arc_num; ++e)
      if (_source[e] != _root && _target[e] != _root) c += _flow[e] * _cost[e];

    return c;
  }

  double totalFlow() const {
    double tot_flow = 0.0;
    for (int e = _dummy_arc; e < _arc_num; ++e)
      if (_source[e] != _root && _target[e] != _root) tot_flow += _flow[e];

    return tot_flow;
  }

  double dummyFlow() const {
    double tot_flow = 0.0;
    for (int e = 0; e < _dummy_arc; ++e)
      if (_source[e] == _root) tot_flow += _flow[e];

    return tot_flow;
  }

  double updateDummyCost(double value) {
    for (int e = 0; e < _dummy_arc; ++e)
      if (_supply[e] >= 0.0)
        _cost[e] = 0.0;
      else
        _cost[e] = value;
  }

  // Recompute potentials
  void recomputePotential() {
    int j = _thread[_root];
    _pi[_root] = 0.0;

    while (j != _root) {
      int e = _pred[j];
      if (j == _target[e])
        _pi[_target[e]] = _pi[_source[e]] + _cost[e];
      else
        _pi[_source[e]] = _pi[_target[e]] - _cost[e];

      j = _thread[j];
    }
  }

  // Potential of node n
  double potential(int n) const { return _pi[n]; }

  // flow on arc e;
  double arcFlow(int e) const { return _flow[e]; }

  // Runtime in milliseconds
  double runtime() const { return _runtime; }

  // Number of iterations of simplex algorithms
  uint64_t iterations() const { return _iterations; }

  // Set basic parameters
  void setTimelimit(double t) {
    _timelimit = t;
    //    PRINT("INFO: change <timelimit> to %f\n", t);
  }
  void setOptTolerance(int o) {
    _opt_tolerance = o;
    //    PRINT("INFO: change <opt_tolerance> to %f\n", o);
  }
  void setVerbosity(std::string v) {
    _verbosity = v;
    if (v == DOT_VAL_DEBUG) N_IT_LOG = 100000;
    if (v == DOT_VAL_INFO) N_IT_LOG = 10000000;
    if (v == DOT_VAL_SILENT) N_IT_LOG = 0;
    //    PRINT("INFO: change <verbosity> to %s\n", v.c_str());
  }

  // Check feasibility
  ProblemType checkFeasibility() {
    for (int e = 0; e != _dummy_arc; ++e)
      if (fabs(_flow[e]) > 1e-09)
        throw std::runtime_error(
            "ERROR 3: flow on dummy arcs: " + std::to_string(_flow[e]) + "\n");

    return ProblemType::OPTIMAL;
  }

  // Reserve memory for arcs in the simplex network
  void reserveArcMemory(size_t s) {
    auto o = _source.size();
    _source.reserve(o + s);
    _target.reserve(o + s);
    _cost.reserve(o + s);

    _flow.reserve(o + s);
    _state.reserve(o + s);
  }

  // Reserve memory for arcs in the simplex network
  void resizeArcMemory(int s) {
    int o = static_cast<int>(_source.size());
    _source.resize(o + s, -1);
    _target.resize(o + s, -1);
    _cost.resize(o + s, -1);

    _flow.resize(o + s, 0);
    _state.resize(o + s, STATE_LOWER);
  }

 private:
  // Initialize internal data structures
  bool init() {
    if (_node_num == 0) return false;

    // Check the sum of supply values
    _sum_supply = 0.0;
    for (int i = 0; i != _node_num; ++i) _sum_supply += _supply[i];

    // Initialize artifical cost
    double ART_COST;
    if (std::numeric_limits<double>::is_exact) {
      ART_COST = (std::numeric_limits<double>::max)() / 2 + 1;
    } else {
      ART_COST = 0.0;
      for (int i = _dummy_arc; i != _arc_num; ++i) {
        if (_cost[i] > ART_COST) ART_COST = _cost[i];
      }
      ART_COST = (ART_COST + 1) * _node_num;
    }

    // Set data for the artificial root node
    _root = _node_num;
    _parent[_root] = -1;
    _pred[_root] = -1;
    _thread[_root] = 0;
    _rev_thread[0] = _root;
    _succ_num[_root] = _node_num + 1;
    _last_succ[_root] = _root - 1;
    _supply[_root] = -_sum_supply;
    _pi[_root] = 0.0;

    for (int u = 0, e = 0; u != _node_num; ++u, ++e) {
      _parent[u] = _root;
      _pred[u] = e;
      _thread[u] = u + 1;
      _rev_thread[u + 1] = u;
      _succ_num[u] = 1;
      _last_succ[u] = u;
      _state[e] = STATE_TREE;
      if (_supply[u] >= 0) {
        _pred_dir[u] = DIR_UP;
        _pi[u] = 0.0;
        _source[e] = u;
        _target[e] = _root;
        _flow[e] = _supply[u];
        _cost[e] = 0.0;
      } else {
        _pred_dir[u] = DIR_DOWN;
        _pi[u] = ART_COST;
        _source[e] = _root;
        _target[e] = u;
        _flow[e] = -_supply[u];
        _cost[e] = ART_COST;
      }
    }

    return true;
  }

  bool pricing() {
    double min = -_opt_tolerance;
    int e = _next_arc;
    int M = _arc_num;
    int M0 = _dummy_arc;
    int e_max = std::min(M, e + _block_size);

    while (e < M) {
      eati_search_dbl(e, e_max, &_state[0], &_cost[0], &_pi[0], &_source[0],
                      &_target[0], &min, &in_arc);

      if (min < std::nextafter(-_opt_tolerance, -DBL_MAX)) {
        _next_arc = in_arc;
        return true;
      }

      e = e_max;
      e_max = std::min(M, e + _block_size);
    }

    e = M0;
    e_max = std::min(_next_arc, e + _block_size);

    while (e < _next_arc) {
      eati_search_dbl(e, e_max, &_state[0], &_cost[0], &_pi[0], &_source[0],
                      &_target[0], &min, &in_arc);

      if (min < std::nextafter(-_opt_tolerance, -DBL_MAX)) {
        _next_arc = in_arc;
        return true;
      }

      e = e_max;
      e_max = std::min(_next_arc, e + _block_size);
    }

    return false;
  }

  // Find the join node
  void findJoinNode() {
    int u = _source[in_arc];
    int v = _target[in_arc];
    while (u != v) {
      if (_succ_num[u] < _succ_num[v]) {
        u = _parent[u];
      } else {
        v = _parent[v];
      }
    }
    join = u;
  }

  // Find the leaving arc of the cycle and returns true if the
  // leaving arc is not the same as the entering arc
  bool findLeavingArc() {
    // Initialize first and second nodes according to the direction
    // of the cycle
    int first, second;
    first = _source[in_arc];
    second = _target[in_arc];

    delta = MAX;
    int result = 0;
    double d;
    int e;

    // Search the cycle form the first node to the join node
    for (int u = first; u != join; u = _parent[u]) {
      e = _pred[u];
      d = _flow[e];
      if (_pred_dir[u] == DIR_DOWN) d = INF - d;

      if (d < delta) {
        delta = d;
        u_out = u;
        result = 1;
      }
    }

    // Search the cycle form the second node to the join node
    for (int u = second; u != join; u = _parent[u]) {
      e = _pred[u];
      d = _flow[e];
      if (_pred_dir[u] == DIR_UP) d = INF - d;

      if (d <= delta) {
        delta = d;
        u_out = u;
        result = 2;
      }
    }

    if (result == 1) {
      u_in = first;
      v_in = second;
    } else {
      u_in = second;
      v_in = first;
    }

    return result != 0;
  }

  // Change _flow and _state vectors
  void changeFlow() {
    // Augment along the cycle
    if (delta > 0.0) {
      _flow[in_arc] += delta;

      //#pragma omp parallel sections num_threads(2)
      {
        //#pragma omp section
        {
          for (int u = _source[in_arc]; u != join; u = _parent[u]) {
            _flow[_pred[u]] -= _pred_dir[u] * delta;
          }
        }
        //#pragma omp section
        {
          for (int u = _target[in_arc]; u != join; u = _parent[u]) {
            _flow[_pred[u]] += _pred_dir[u] * delta;
          }
        }
      }
    }

    // Update the state of the entering and leaving arcs
    _state[in_arc] = STATE_TREE;
    _state[_pred[u_out]] = STATE_LOWER;
  }

  // Update the tree structure
  void updateTreeStructure() {
    int old_rev_thread = _rev_thread[u_out];
    int old_succ_num = _succ_num[u_out];
    int old_last_succ = _last_succ[u_out];
    v_out = _parent[u_out];

    // Check if u_in and u_out coincide
    if (u_in == u_out) {
      // Update _parent, _pred, _pred_dir
      _parent[u_in] = v_in;
      _pred[u_in] = in_arc;
      _pred_dir[u_in] = u_in == _source[in_arc] ? DIR_UP : DIR_DOWN;

      // Update _thread and _rev_thread
      if (_thread[v_in] != u_out) {
        int after = _thread[old_last_succ];
        _thread[old_rev_thread] = after;
        _rev_thread[after] = old_rev_thread;
        after = _thread[v_in];
        _thread[v_in] = u_out;
        _rev_thread[u_out] = v_in;
        _thread[old_last_succ] = after;
        _rev_thread[after] = old_last_succ;
      }
    } else {
      // Handle the case when old_rev_thread equals to v_in
      // (it also means that join and v_out coincide)
      int thread_continue =
          old_rev_thread == v_in ? _thread[old_last_succ] : _thread[v_in];

      // Update _thread and _parent along the stem nodes (i.e. the nodes
      // between u_in and u_out, whose parent have to be changed)
      int stem = u_in;              // the current stem node
      int par_stem = v_in;          // the new parent of stem
      int next_stem;                // the next stem node
      int last = _last_succ[u_in];  // the last successor of stem
      int before, after = _thread[last];
      _thread[v_in] = u_in;
      _dirty_revs.clear();
      _dirty_revs.push_back(v_in);
      while (stem != u_out) {
        // Insert the next stem node into the thread list
        next_stem = _parent[stem];
        _thread[last] = next_stem;
        _dirty_revs.push_back(last);

        // Remove the subtree of stem from the thread list
        before = _rev_thread[stem];
        _thread[before] = after;
        _rev_thread[after] = before;

        // Change the parent node and shift stem nodes
        _parent[stem] = par_stem;
        par_stem = stem;
        stem = next_stem;

        // Update last and after
        last = _last_succ[stem] == _last_succ[par_stem] ? _rev_thread[par_stem]
                                                        : _last_succ[stem];
        after = _thread[last];
      }
      _parent[u_out] = par_stem;
      _thread[last] = thread_continue;
      _rev_thread[thread_continue] = last;
      _last_succ[u_out] = last;

      // Remove the subtree of u_out from the thread list except for
      // the case when old_rev_thread equals to v_in
      if (old_rev_thread != v_in) {
        _thread[old_rev_thread] = after;
        _rev_thread[after] = old_rev_thread;
      }

      // Update _rev_thread using the new _thread values
      for (int i = 0; i != int(_dirty_revs.size()); ++i) {
        int u = _dirty_revs[i];
        _rev_thread[_thread[u]] = u;
      }

      // Update _pred, _pred_dir, _last_succ and _succ_num for the
      // stem nodes from u_out to u_in
      int tmp_sc = 0, tmp_ls = _last_succ[u_out];
      for (int u = u_out, p = _parent[u]; u != u_in; u = p, p = _parent[u]) {
        _pred[u] = _pred[p];
        _pred_dir[u] = -_pred_dir[p];
        tmp_sc += _succ_num[u] - _succ_num[p];
        _succ_num[u] = tmp_sc;
        _last_succ[p] = tmp_ls;
      }
      _pred[u_in] = in_arc;
      _pred_dir[u_in] = u_in == _source[in_arc] ? DIR_UP : DIR_DOWN;
      _succ_num[u_in] = old_succ_num;
    }

    // Update _last_succ from v_in towards the root
    int up_limit_out = _last_succ[join] == v_in ? join : -1;
    int last_succ_out = _last_succ[u_out];
    for (int u = v_in; u != -1 && _last_succ[u] == v_in; u = _parent[u]) {
      _last_succ[u] = last_succ_out;
    }

    // Update _last_succ from v_out towards the root
    if (join != old_rev_thread && v_in != old_rev_thread) {
      for (int u = v_out; u != up_limit_out && _last_succ[u] == old_last_succ;
           u = _parent[u]) {
        _last_succ[u] = old_rev_thread;
      }
    } else if (last_succ_out != old_last_succ) {
      for (int u = v_out; u != up_limit_out && _last_succ[u] == old_last_succ;
           u = _parent[u]) {
        _last_succ[u] = last_succ_out;
      }
    }

    // Update _succ_num from v_in to join
    for (int u = v_in; u != join; u = _parent[u]) {
      _succ_num[u] += old_succ_num;
    }
    // Update _succ_num from v_out to join
    for (int u = v_out; u != join; u = _parent[u]) {
      _succ_num[u] -= old_succ_num;
    }
  }

  // Update potentials in the subtree that has been moved
  void updatePotential() {
    double sigma = _pi[v_in] - _pi[u_in] - _pred_dir[u_in] * _cost[in_arc];
    int end = _thread[_last_succ[u_in]];

    if (sigma != 0)
      for (int u = u_in; u != end; u = _thread[u]) _pi[u] += sigma;
  }

  // Execute the algorithm
  ProblemType start(PivotRule pivot_rule) {
    // Select the pivot rule implementation
    switch (pivot_rule) {
      case PivotRule::BLOCK_SEARCH:
        return start<BlockSearchPivotRule>();
    }
    return ProblemType::INFEASIBLE;  // avoid warning
  }

  template <typename PivotRuleImpl>
  ProblemType start() {
    auto start_tt = std::chrono::steady_clock::now();
    PivotRuleImpl pivot(*this);

    // Benchmarking

    // Execute the Network Simplex algorithm
    t1 = 0.0;
    while (true) {
      _iterations++;

      bool stop = pricing();
      // bool stop = pivot.findEnteringArc();
      fprintf(stdout, "beep1");

      if (!stop) break;

      findJoinNode();
      findLeavingArc();
      changeFlow();
      updateTreeStructure();
      updatePotential();

      fprintf(stdout, "beep");
    }

    auto end_t = std::chrono::steady_clock::now();
    _runtime += double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_tt)
                           .count()) /
                1000;

    if (_verbosity == DOT_VAL_DEBUG)
      PRINT(
          "NetSIMPLEX outer loop | enter: %.3f, join: %.3f, leave: %.3f, "
          "change: %.3f, tree: %.3f, "
          "potential: %.3f, runtime: %.3f\n",
          t1, t2, t3, t4, t5, t6, _runtime);

    return ProblemType::OPTIMAL;
  }
};
}  // namespace DOT
