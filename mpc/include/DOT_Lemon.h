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

// Main network simplex implementation with double for cost, capacity, and
// supplies
class NetSimplexCapacity {
 public:
  /// The type of the flow amounts, capacity bounds and supply values
  typedef double Value;
  /// The type of the arc costs
  typedef double Cost;

 private:
  typedef std::vector<int> IntVector;
  typedef std::vector<Value> ValueVector;
  typedef std::vector<Cost> CostVector;
  typedef std::vector<signed char> CharVector;

  // State constants for arcs
  const int STATE_UPPER = -1;
  const int STATE_TREE = 0;
  const int STATE_LOWER = 1;

  // Direction constants for tree arcs
  const int DIR_DOWN = -1;
  const int DIR_UP = 1;

  // Data related to the underlying digraph
  int _node_num;
  int _arc_num;

  int _fixed_arc;
  int _dummy_arc;  // Arc id where begin the basic arcs
  int _next_arc;

  int _all_arc_num;
  int _search_arc_num;

  // Parameters of the problem
  Value _sum_supply;

  // Data structures for storing the digraph
  IntVector _source;
  IntVector _target;

  // Node and arc data
  ValueVector _upper;
  ValueVector _cap;

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
  CharVector _state;
  IntVector _dirty_revs;

  int _root;

  // Temporary data used in the current pivot iteration
  int in_arc, join, u_in, v_in, u_out, v_out;
  Value delta;

  const Value MAX;
  const Value INF;

  double _runtime;

  double _timelimit;
  double _opt_tolerance;

  size_t _iterations;

 public:
  // Implementation of the Block Search pivot rule
  class BlockSearchPivotRule {
   private:
    // References to the NetworkSimplex class
    const IntVector& _source;
    const IntVector& _target;
    const CostVector& _cost;
    const CharVector& _state;
    const CostVector& _pi;
    int& _in_arc;
    int _arc_num;
    int _dummy_arc;

    // Pivot rule data
    int _block_size;
    int _next_arc;

    // Negative eps
    const double negeps;

   public:
    // Constructor
    BlockSearchPivotRule(NetSimplexCapacity& ns)
        : _source(ns._source),
          _target(ns._target),
          _cost(ns._cost),
          _state(ns._state),
          _pi(ns._pi),
          _in_arc(ns.in_arc),
          _arc_num(ns._arc_num),
          _dummy_arc(ns._dummy_arc),
          _next_arc(ns._next_arc),
          negeps(std::nextafter(-ns._opt_tolerance, -0.0)) {
      // The main parameters of the pivot rule
      const double BLOCK_SIZE_FACTOR = 1.0;
      const int MIN_BLOCK_SIZE =
          20;  // THIS VALUE IS IMPORTANT AND IT SHOULD BE A PARAMETER (!)

      _block_size =
          std::max(int(BLOCK_SIZE_FACTOR * std::sqrt(double(_dummy_arc))),
                   MIN_BLOCK_SIZE);
    }

    // Find next entering arc
    bool findEnteringArc() {
      Cost c, min = negeps;
      int cnt = _block_size;
      int e;
      for (e = _next_arc; e != _arc_num; ++e) {
        c = _state[e] * (_cost[e] + _pi[_source[e]] - _pi[_target[e]]);
        if (c < min) {
          min = c;
          _in_arc = e;
        }
        if (--cnt == 0) {
          if (min < negeps) goto search_end;
          cnt = _block_size;
        }
      }
      for (e = _dummy_arc; e != _next_arc; ++e) {
        c = _state[e] * (_cost[e] + _pi[_source[e]] - _pi[_target[e]]);
        if (c < min) {
          min = c;
          _in_arc = e;
        }
        if (--cnt == 0) {
          if (min < negeps) goto search_end;
          cnt = _block_size;
        }
      }
      if (min >= negeps) return false;

    search_end:
      _next_arc = _in_arc;
      return true;
    }

  };  // class BlockSearchPivotRule

 public:
  NetSimplexCapacity(const char INIT, int node_num, int arc_num)
      : _node_num(node_num),
        _arc_num(0),
        _root(-1),
        in_arc(-1),
        join(-1),
        u_in(-1),
        v_in(-1),
        u_out(-1),
        v_out(-1),
        MAX((std::numeric_limits<Value>::max)()),
        INF(std::numeric_limits<Value>::has_infinity
                ? std::numeric_limits<Value>::infinity()
                : MAX),
        _runtime(0.0) {
    // Check the number types
    if (!std::numeric_limits<Value>::is_signed)
      throw std::runtime_error(
          "The flow type of NetworkSimplex must be signed");
    if (!std::numeric_limits<Cost>::is_signed)
      throw std::runtime_error(
          "The cost type of NetworkSimplex must be signed");

    // Reset data structures
    int all_node_num = _node_num + 1;
    _supply.resize(all_node_num, 0);
    _pi.resize(all_node_num);

    _parent.resize(all_node_num);
    _pred.resize(all_node_num);
    _pred_dir.resize(all_node_num);
    _thread.resize(all_node_num);
    _rev_thread.resize(all_node_num);
    _succ_num.resize(all_node_num);
    _last_succ.resize(all_node_num);

    // 2*n arcs from nodes to root and from root to node;
    // 2*n-1 nodes in a basic solution
    int max_arc_num = 0;
    if (INIT == 'F')  // Full
      max_arc_num = 2 * _node_num + arc_num + 1;

    if (INIT == 'E')  // Empty, for Column Generation
      max_arc_num = 4 * _node_num + 1;

    _source.reserve(max_arc_num);
    _target.reserve(max_arc_num);

    _upper.reserve(max_arc_num);
    _cap.reserve(max_arc_num);

    _cost.reserve(max_arc_num);
    _flow.reserve(max_arc_num);
    _state.reserve(max_arc_num);

    _source.resize(_node_num);
    _target.resize(_node_num);

    _cost.resize(_node_num, 0);
    _flow.resize(_node_num, 0);
    _state.resize(_node_num, STATE_LOWER);

    _upper.resize(_node_num, MAX);
    _cap.resize(_node_num, MAX);

    // Dummy arcs for eavery node to root node
    _fixed_arc = 0;
    _dummy_arc = _node_num;
    _arc_num = _node_num;
    _next_arc = _dummy_arc;

    // Interal parameters
    _timelimit = std::numeric_limits<double>::max();
    _opt_tolerance = 1e-06;
  }

  ProblemType run() {
    _runtime = 0.0;
    _iterations = 0;
    // Reset arc variables
    for (int e = 0; e < _arc_num; ++e) {
      _state[e] = STATE_LOWER;
      _flow[e] = 0.0;
    }

    if (!init()) return ProblemType::INFEASIBLE;
    return start();
  }

  ProblemType reRun() { return start(); }

  int num_arcs() const { return int(_source.size()) - _dummy_arc; }

  int num_nodes() const { return _node_num; }

  void addNode(int i, Value b) { _supply[i] = b; }

  size_t addArc(int a, int b, Cost c, Value u) {
    size_t idx = _source.size();
    _source.push_back(a);
    _target.push_back(b);
    _cost.push_back(c);

    // Recheck
    _cap.push_back(u);
    _upper.push_back(u);

    _flow.push_back(0);
    _state.push_back(STATE_LOWER);

    _arc_num++;
    return idx;
  }

  // Change the cost to a single arc
  void setArcCost(size_t idx, Cost value) { _cost[idx] = value; }

  void fixArcs(void) { _fixed_arc = _arc_num; }

  void setArc(size_t idx, int a, int b, Cost c, Value u) {
    _source[_dummy_arc + idx] = a;
    _target[_dummy_arc + idx] = b;
    _cost[_dummy_arc + idx] = c;

    _cap[_dummy_arc + idx] = u;
    _upper[_dummy_arc + idx] = u;

    _flow[_dummy_arc + idx] = 0;
    _state[_dummy_arc + idx] = STATE_LOWER;

    _arc_num++;
  }

  // MANCA UPDATE ARCS
  int updateArcs(const Vars& as) {
    int new_arc = 0;
    size_t idx = 0;
    size_t idx_max = as.size();

    int e = _fixed_arc;
    int e_max = _arc_num;

    // Store the new arc variable, replacing an used arc,
    // out of the basis and with positive reduced cost
    for (; idx < idx_max; ++idx) {
      while (e < e_max) {
        // Replace useless variables with new variables
        if (_state[e] == STATE_LOWER &&
            (_cost[e] + _pi[_source[e]] - _pi[_target[e]] > 1e-09))
          break;
        ++e;
      }
      if (e >= e_max) break;
      _source[e] = as[idx].a;
      _target[e] = as[idx].b;
      _cost[e] = as[idx].c;
      _cap[e] = MAX;
      _upper[e] = MAX;
      if (new_arc == 0) _next_arc = e;
      new_arc++;
    }

    for (; idx < idx_max; ++idx) {
      addArc(as[idx].a, as[idx].b, as[idx].c, MAX);
      if (new_arc == 0) _next_arc = e;
      new_arc++;
    }

    return new_arc;
  }

  Cost totalCost() const {
    Cost c = 0;
    for (int e = _dummy_arc; e < _arc_num; ++e)
      if (_source[e] != _root && _target[e] != _root) {
        c += _flow[e] * _cost[e];
      }

    return c;
  }

  Cost totalFlow() const {
    Cost tot_flow = 0;
    for (int e = _dummy_arc; e < _arc_num; ++e)
      if (_source[e] != _root && _target[e] != _root) tot_flow += _flow[e];

    return tot_flow;
  }

  // Potential of node n
  Cost potential(int n) const { return _pi[n]; }

  // Runtime in milliseconds
  double runtime() const { return _runtime; }

  // Number of iterations of simplex algorithms
  int iterations() const { return _iterations; }

  // Set basic parameters
  void setTimelimit(double t) { _timelimit = t; }

  void setOptTolerance(double o) { _opt_tolerance = o; }

  // Check feasibility
  ProblemType checkFeasibility() {
    for (int e = 0; e != _dummy_arc; ++e)
      if (fabs(_flow[e]) > 1e-09 && _flow[e] <= _cap[e])
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

    _cap.reserve(o + s);
    _upper.reserve(o + s);

    _flow.reserve(o + s);
    _state.reserve(o + s);
  }

  // Reserve memory for arcs in the simplex network
  void resizeArcMemory(int s) {
    int o = static_cast<int>(_source.size());
    _source.resize(o + s, -1);
    _target.resize(o + s, -1);
    _cost.resize(o + s, -1);

    _cap.resize(o + s);
    _upper.resize(o + s);

    _flow.resize(o + s, 0);
    _state.resize(o + s, STATE_LOWER);
  }

 private:
  // Initialize internal data structures
  bool init() {
    if (_node_num == 0) return false;

    // Check the sum of supply values
    _sum_supply = 0;
    for (int i = 0; i != _node_num; ++i) _sum_supply += _supply[i];

    if (_sum_supply > 0) return false;

    // Initialize artifical cost
    Cost ART_COST;
    if (std::numeric_limits<Cost>::is_exact) {
      ART_COST = (std::numeric_limits<Cost>::max)() / 2 + 1;
    } else {
      ART_COST = 0;
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
    _pi[_root] = 0;

    // Add artificial arcs and initialize the spanning tree data structure
    if (_sum_supply == 0) {
      // EQ supply constraints
      for (int u = 0, e = 0; u != _node_num; ++u, ++e) {
        _parent[u] = _root;
        _pred[u] = e;
        _thread[u] = u + 1;
        _rev_thread[u + 1] = u;
        _succ_num[u] = 1;
        _last_succ[u] = u;
        _cap[e] = INF;
        _state[e] = STATE_TREE;
        if (_supply[u] >= 0) {
          _pred_dir[u] = DIR_UP;
          _pi[u] = 0;
          _source[e] = u;
          _target[e] = _root;
          _flow[e] = _supply[u];
          _cost[e] = 0;
        } else {
          _pred_dir[u] = DIR_DOWN;
          _pi[u] = ART_COST;
          _source[e] = _root;
          _target[e] = u;
          _flow[e] = -_supply[u];
          _cost[e] = ART_COST;
        }
      }
    } else {
      // GEQ supply constraints
      int f = _arc_num + _node_num;
      for (int u = 0, e = _arc_num; u != _node_num; ++u, ++e) {
        _parent[u] = _root;
        _thread[u] = u + 1;
        _rev_thread[u + 1] = u;
        _succ_num[u] = 1;
        _last_succ[u] = u;
        if (_supply[u] <= 0) {
          _pred_dir[u] = DIR_DOWN;
          _pi[u] = 0;
          _pred[u] = e;
          _source[e] = _root;
          _target[e] = u;
          _cap[e] = INF;
          _flow[e] = -_supply[u];
          _cost[e] = 0;
          _state[e] = STATE_TREE;
        } else {
          _pred_dir[u] = DIR_UP;
          _pi[u] = -ART_COST;
          _pred[u] = f;
          _source[f] = u;
          _target[f] = _root;
          _cap[f] = INF;
          _flow[f] = _supply[u];
          _state[f] = STATE_TREE;
          _cost[f] = ART_COST;
          _source[e] = _root;
          _target[e] = u;
          _cap[e] = INF;
          _flow[e] = 0;
          _cost[e] = 0;
          _state[e] = STATE_LOWER;
          ++f;
        }
      }
      _all_arc_num = f;
    }

    return true;
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
    if (_state[in_arc] == STATE_LOWER) {
      first = _source[in_arc];
      second = _target[in_arc];
    } else {
      first = _target[in_arc];
      second = _source[in_arc];
    }
    delta = _cap[in_arc];
    int result = 0;
    Value c, d;
    int e;

    // Search the cycle form the first node to the join node
    for (int u = first; u != join; u = _parent[u]) {
      e = _pred[u];
      d = _flow[e];
      if (_pred_dir[u] == DIR_DOWN) {
        c = _cap[e];
        d = c >= MAX ? INF : c - d;
      }
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
      if (_pred_dir[u] == DIR_UP) {
        c = _cap[e];
        d = c >= MAX ? INF : c - d;
      }
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
  void changeFlow(bool change) {
    // Augment along the cycle
    if (delta > 0) {
      Value val = _state[in_arc] * delta;
      _flow[in_arc] += val;
      for (int u = _source[in_arc]; u != join; u = _parent[u])
        _flow[_pred[u]] -= _pred_dir[u] * val;
      for (int u = _target[in_arc]; u != join; u = _parent[u])
        _flow[_pred[u]] += _pred_dir[u] * val;

      // Update the state of the entering and leaving arcs
      if (change) {
        _state[in_arc] = STATE_TREE;
        _state[_pred[u_out]] =
            (_flow[_pred[u_out]] == 0) ? STATE_LOWER : STATE_UPPER;
      } else {
        _state[in_arc] = -_state[in_arc];
      }
    }
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
    Cost sigma = _pi[v_in] - _pi[u_in] - _pred_dir[u_in] * _cost[in_arc];
    int end = _thread[_last_succ[u_in]];
    for (int u = u_in; u != end; u = _thread[u]) {
      _pi[u] += sigma;
    }
  }

  ProblemType start(void) {
    auto start_tt = std::chrono::steady_clock::now();

    BlockSearchPivotRule pivot(*this);

    // Execute the Network Simplex algorithm
    while (true) {
      bool stop = pivot.findEnteringArc();
      if (!stop) break;

      findJoinNode();
      bool change = findLeavingArc();
      if (delta >= MAX) return ProblemType::UNBOUNDED;
      changeFlow(change);
      if (change) {
        updateTreeStructure();
        updatePotential();
      }

      // Add as log file
      _iterations++;
      if (_iterations % 1000 == 0)
        fprintf(stdout, "%d - %.4f\n", _iterations, totalCost());
    }

    auto end_t = std::chrono::steady_clock::now();
    _runtime += double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_tt)
                           .count()) /
                1000;

    return ProblemType::OPTIMAL;
  }
};

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
    capLB = &data[_m];
    capUB = &data[2 * _m];
    supply = &data[3 * m];
  }

  // Reset memory for reusing the allocated memory
  void reset(void) {
    memset(arcs, 0, size_arcs);
    memset(data, 0, size_data);
  }

  // Add node
  void addNode(size_t i, double b) { supply[i] = b; }

  // Add arc
  void addArc(size_t i, size_t j, double _cost, double _capLB, double _capUB) {
    head[idx] = i;
    tail[idx] = j;

    cost[idx] = _cost;
    capLB[idx] = _capLB;
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
          n = std::stol(buf);
          buf = strchr(buf, sep);
          ++buf;
          m = std::stol(buf);
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

    while (buf != NULL) {
      if (buf[0] == 'n') {
        buf = strchr(buf, sep);
        ++buf;

        size_t _idx = std::stol(buf);

        buf = strchr(buf, sep);
        ++buf;
        double _s = std::stod(buf);

        buf = strchr(buf, eol);
        if (buf == NULL) break;
        ++buf;

        addNode(_idx - 1, _s);
      } else {
        if (buf[0] == 'a') {
          buf = strchr(buf, sep);
          ++buf;

          size_t i = std::stol(buf);

          buf = strchr(buf, sep);
          ++buf;
          size_t j = std::stol(buf);

          buf = strchr(buf, sep);
          ++buf;
          double _lb = std::stod(buf);

          buf = strchr(buf, sep);
          ++buf;
          double _ub = std::stod(buf);

          buf = strchr(buf, sep);
          ++buf;
          double _c = std::stod(buf);
          addArc(i - 1, j - 1, _c, _lb, _ub);
          if (idx == m) return;

          buf = strchr(buf, eol);
          ++buf;
        }
      }
    }
  }

  // Solve the instance
  double solve() {
    NetSimplexCapacity simplex('F', n, m);
    simplex.setTimelimit(3600);
    simplex.setOptTolerance(1.0);

    for (size_t i = 0; i < n; i++) simplex.addNode((int)i, supply[i]);

    for (size_t i = 0; i < m; i++)
      simplex.addArc((int)head[i], (int)tail[i], cost[i], capUB[i]);

    simplex.run();

    double r = simplex.totalCost();

    size_t api_it = simplex.iterations();
    double api_time = simplex.runtime();

    fprintf(stdout, "CHECK SOL %.3f ==> %d, it: %d, time: %.3f\n", n, r, api_it,
            api_time);

    return r;
  }

  // Dump the graph to a file
  void dump(void) const {
    for (size_t u = 0; u < n; u++) fprintf(stdout, "%ld %f\n", u, supply[u]);

    for (size_t e = 0; e < m; e++)
      fprintf(stdout,
              "(head: %ld, tail: %ld, lb: %.3f, ub: %.3f, cost: %.3f)\n",
              head[e], tail[e], capLB[e], capUB[e], cost[e]);
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

  // Keep contigous memory
  size_t size_data;
  double* data;
  double* cost;    // one for arcs
  double* capLB;   // one for arcs
  double* capUB;   // one for arcs
  double* supply;  // one for node
};
}  // namespace DOT
