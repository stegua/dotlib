/**
 * @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include "DOT_Commons.h"
#include <vector>
using std::vector;

#include <algorithm>
#include <numeric>

namespace DOT {

class EapiSimplex {
public:
  // State constants for arcs
  const bool STATE_TREE = false;
  const bool STATE_LOWER = true;

  // Direction constants for tree arcs
  const char DIR_DOWN = -1;
  const char DIR_UP = 1;

  // Parameters of the problem
  int64_t total_mass;

  // Arc data
  vector<int> source;
  vector<int> target;
  vector<int> flow;
  vector<int> cost;

  // Node data
  vector<int> balance;
  vector<int> pi;

  // Data for storing the spanning tree basis
  vector<int> father;
  vector<int> fatherArc;
  vector<int> son;
  vector<int> brother;
  vector<int> sister;
  vector<int> successors;

  vector<char> fatherArcDir;
  vector<bool> state;

  int root;

  // Data related to the underlying digraph
  int N;
  int M;
  int M0;        // 2*N dummy arcs
  int tot_arcs;  // 2*N + M arcs (2*N dummy arcs)
  int tot_nodes; // N + 1 root

  int in_arc;    // Arc entering the basis
  int u_in;      // Source of entering node
  int v_in;      // Target of entering node
  int join_node; // Node closing the cycle
  int delta;     // amount of flow change in the cycle

public:
  // Standard c'tor
  EapiSimplex(int _n, int _m)
      : root(_n), N(_n), M(_m), M0(2 * _n), tot_arcs(2 * N + M),
        tot_nodes(N + 1) {
    // Arcs
    source.reserve(tot_arcs);
    target.reserve(tot_arcs);
    flow.reserve(tot_arcs);
    cost.reserve(tot_arcs);
    state.reserve(tot_arcs);

    source.resize(M0, 0);
    target.resize(M0, 0);
    cost.resize(M0, 0);
    flow.resize(tot_arcs, 0);
    state.resize(tot_arcs, STATE_LOWER);

    // Node data
    balance.resize(tot_nodes, 0);
    pi.resize(tot_nodes, 0);

    // Data for storing the spanning tree basis
    father.resize(tot_nodes, 0);
    fatherArc.resize(tot_nodes, 0);
    fatherArcDir.resize(tot_nodes, 0);

    son.resize(tot_nodes, 0);
    brother.resize(tot_nodes, 0);
    sister.resize(tot_nodes, 0);
    successors.resize(tot_nodes, 0);
  }

  // Add a node balance
  void addNode(int i, int b) { balance[i] = b; }

  // Add an arc
  int addArc(int i, int j, int c) {
    source.push_back(i);
    target.push_back(j);
    cost.push_back(c);
  }

private:
  // Initialize internal data structures
  bool init() {
    // Start with the basis tree
    int tot_balance = std::accumulate(balance.begin(), balance.end(), 0);

    // Initialize artifical cost
    int ART_COST = 1 + *(std::max_element(cost.begin(), cost.end()));

    // Set data for the artificial root node
    root = N;
    balance[root] = -tot_balance;
    pi[root] = 0;

    father[root] = -1;
    fatherArc[root] = -1;
    brother[root] = -1;
    sister[root] = -1;
    son[root] = 0;
    successors[root] = N;

    // Init spanning tree
    for (int v = 0, e = 0; v < N; ++v, ++e) {
      father[v] = root;
      fatherArc[v] = e;
      fatherArcDir[v] = -1;
      brother[v] = v + 1;
      sister[v] = v - 1;
      son[v] = 0;
      successors[v] = 0;
      state[e] = STATE_TREE;

      if (balance[v] >= 0) {
        fatherArcDir[e] = DIR_UP;
        source[e] = v;
        target[e] = root;
        flow[e] = balance[v];
      } else {
        fatherArcDir[e] = DIR_DOWN;
        pi[v] = ART_COST;
        source[e] = root;
        target[e] = v;
        flow[e] = -balance[v];
        cost[e] = ART_COST;
      }
    }
  }

  // Pricing for negative reduced cost
  bool pricing() {
    int best = 0;

    for (int e = M0; e < M; ++e) {
      // Ma serve lo stato? Se è in base dovrebbe sempre essere positivo!
      int c = state[e] * (cost[e] + pi[source[e]] - pi[target[e]]);
      if (c < best) {
        best = c;
        in_arc = e;
      }
    }

    return best < 0;
  }

  // Look for the arc leaving the basis
  void findJoinNode() {
    int i = source[in_arc];
    int j = target[in_arc];
    int p = i;
    int k = j;
    if (successors[i] < successors[j]) {
      p = j;
      k = i;
    }
    while (successors[p] > successors[k])
      p = father[p];
    while (p != k) {
      p = father[p];
      k = father[k];
    }
    join_node = p;
  }

  // Look for leaving variable
  void findLeavingArc() {
    int first = source[in_arc];
    int second = target[in_arc];

    int MAX = std::numeric_limits<int>::max();
    delta = MAX;
    int result = 0;
    int d;
    int e;

    // Search the cycle form the first node to the join node
    for (int u = first; u != join_node; u = father[u]) {
      e = fatherArc[u];
      d = flow[e];
      if (fatherArcDir[u] == DIR_DOWN)
        d = MAX - d;

      if (d < delta) {
        delta = d;
        join_node = u;
        result = 1;
      }
    }

    // Search the cycle form the second node to the join node
    for (int u = second; u != join_node; u = father[u]) {
      e = fatherArc[u];
      d = flow[e];
      if (fatherArcDir[u] == DIR_UP)
        d = MAX - d;

      if (d <= delta) {
        delta = d;
        join_node = u;
        result = 2;
      }
    }

    // TODO: Check if this is correct?
    if (result == 1) {
      u_in = first;
      v_in = second;
    } else {
      u_in = second;
      v_in = first;
    }
  }

  // Run the Netowrk Simplex algorithm
  int solve() {
    // Init the internal data structures
    init();

    // Loop for negative reduced cost cycles
    while (true) {
      // Look for negative reduced cost variables out of the basis
      bool stop = pricing();
      if (stop)
        break;

      // Look for exiting variable
      findJoinNode();
      findLeavingArc();
    }
  }
};
} // namespace DOT
