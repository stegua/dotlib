/**
 * @fileoverview Copyright (c) 2019-2022, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <vector>

#include "DOT_Commons.h"
using std::vector;

#include <algorithm>
#include <cassert>
#include <numeric>

namespace DOT {

// Change flow along the upper path from node start (u_in or v_in) to stop
// (i.e., the join_node)
inline static void para_flow(int start, int stop, const int64_t *father,
                             int64_t *flow_pred, const int64_t *dir_pred,
                             int64_t delta) {
  for (int u = start; u != stop; u = father[u])
    flow_pred[u] -= dir_pred[u] * delta;
}

static void para3_search(int e, int e_max, const int32_t *state,
                         const int32_t *cost, const int32_t *pi,
                         const int *source, const int *target, int *min,
                         int *in_arc) {
  min[0] = 0;
  int values[4] = {0, 0, 0, 0};

  for (size_t i = e; i < e_max; i += 4) {
    values[0] = state[i] * (cost[i] + pi[source[i]] - pi[target[i]]);
    values[1] =
        state[i + 1] * (cost[i + 1] + pi[source[i + 1]] - pi[target[i + 1]]);
    values[2] =
        state[i + 2] * (cost[i + 2] + pi[source[i + 2]] - pi[target[i + 2]]);
    values[3] =
        state[i + 3] * (cost[i + 3] + pi[source[i + 3]] - pi[target[i + 3]]);

    if (values[0] < min[0]) {
      min[0] = values[0];
      in_arc[0] = i;
    }
    if (values[1] < min[0]) {
      min[0] = values[1];
      in_arc[0] = i + 1;
    }
    if (values[2] < min[0]) {
      min[0] = values[2];
      in_arc[0] = i + 2;
    }
    if (values[3] < min[0]) {
      min[0] = values[3];
      in_arc[0] = i + 3;
    }
  }
}

//// Perform a search for the pricing
// static void para22_search(int e, int e_max, const int32_t *state,
//                          const int32_t *cost, const int32_t *pi,
//                          const int *source, const int *target, int *min,
//                          int *in_arc) {
//  min[0] = 0;
//  const __m128i increment = _mm_set1_epi32(4);
//  __m128i indices = _mm_setr_epi32(e, e + 1, e + 2, e + 3);
//  __m128i minindices = indices;
//  __m128i minvalues = _mm_setr_epi32(0, 0, 0, 0);
//
//  if ((e_max - e) % 4 != 0)
//    exit(-1);
//
//  for (size_t i = e; i < e_max; i += 4) {
//    const __m128i values = _mm_setr_epi32(
//        state[i] * (cost[i] + pi[source[i]] - pi[target[i]]),
//        state[i + 1] * (cost[i + 1] + pi[source[i + 1]] - pi[target[i + 1]]),
//        state[i + 2] * (cost[i + 2] + pi[source[i + 2]] - pi[target[i + 2]]),
//        state[i + 3] * (cost[i + 3] + pi[source[i + 3]] - pi[target[i + 3]]));
//
//    const __m128i lt = _mm_cmplt_epi32(values, minvalues);
//    minindices = _mm_blendv_epi8(minindices, indices, lt);
//    minvalues = _mm_min_epi32(values, minvalues);
//
//    indices = _mm_add_epi32(indices, increment);
//  }
//
//  // find min index in vector result (in an extremely naive way)
//  int32_t values_array[4];
//  uint32_t indices_array[4];
//
//  _mm_storeu_si128((__m128i *)values_array, minvalues);
//  _mm_storeu_si128((__m128i *)indices_array, minindices);
//
//  for (int i = 0; i < 4; i++) {
//    if (values_array[i] < min[0]) {
//      min[0] = values_array[i];
//      in_arc[0] = (int)indices_array[i];
//    }
//  }
//}

static void para_search(int e, int e_max, const int *state, const int *cost,
                        const int *pi, const int *source, const int *target,
                        int *min, int *in_arc) {
  // min[0] = 0;
  for (int i = e; i < e_max; ++i) {
    int c = state[i] * (cost[i] + pi[source[i]] - pi[target[i]]);
    if (c < min[0]) {
      min[0] = c;
      in_arc[0] = i;
    }
  }
}

// Change flow along the upper path from node start (u_in or v_in) to stop
// (i.e., the join_node)
static void para_flow(int start, int stop, const int *father, int *flow_pred,
                      const int *dir_pred, int delta) {
  for (int u = start; u != stop; u = father[u])
    flow_pred[u] -= dir_pred[u] * delta;
}

// Main class
class EapiSimplex {
 public:
  // State constants for arcs
  const int STATE_TREE = 0;
  const int STATE_LOWER = 1;

  // Direction constants for tree arcs
  const int DIR_DOWN = -1;
  const int DIR_UP = 1;

  // Pricing startegy
  const int BLOCK_SIZE_FACTOR = 4;
  const int MIN_BLOCK_SIZE = 20;
  int _block_size;
  int _next_arc;

  // Parameters of the problem
  int64_t total_mass;

  // Node data
  vector<int> balance;

  // Data for storing the spanning tree basis
  vector<int> data;
  int *father, *son, *brother, *successors, *pi;

  // Only flow and direction for predecessor
  int *flow_pred, *edge_pred, *dir_pred;

  // Arc data
  int *source, *target, *cost, *state;

  //// Only information on arc

  int root;

  // Data related to the underlying digraph
  int N;
  int M;
  int M0;         // 2*N dummy arcs
  int tot_nodes;  // N + 1 root

  int in_arc;     // Arc entering the basis
  int out_arc;    // Arc leaving the basis
  int join_node;  // Node closing the cycle
  int delta;      // amount of flow change in the cycle

  int q;           // Node with entering arc is the leaving arc
  int u_in, v_in;  // First and second node of the entering arc

  double _runtime;
  double _time_pricing;
  double _time_update_duals;
  double _time_update_basis;
  int _iterations;

  // Smallest tree to update the dual multipliers
  bool small_right;

  // Thread pool
  // https://github.com/bshoshany/thread-pool
  //  thread_pool thpool;
  // Parallele pricing
  const int np;
  /*vector<int> arcmin;
  vector<int> arcin;*/

 public:
  // Standard c'tor
  EapiSimplex(int _n, int _m, int _np = 8)
      : root(_n),
        N(_n),
        M(_n),
        M0(_n),
        tot_nodes(N + 1),
        q(-1),
        _runtime(0.0),
        _time_pricing(0),
        _iterations(0),
        np(_np) {
    // Node data
    balance.resize(tot_nodes, 0);

    // Data for storing the spanning tree basis
    data.resize(tot_nodes * 8 + (M + _m) * 4, 0);
    father = &data[0];
    son = &data[tot_nodes];
    brother = &data[tot_nodes * 2];
    successors = &data[tot_nodes * 3];
    pi = &data[tot_nodes * 4];

    dir_pred = &data[tot_nodes * 5];
    flow_pred = &data[tot_nodes * 6];
    edge_pred = &data[tot_nodes * 7];

    // Arcs
    int M4 = M + _m;
    // 4 * (int((M + _m) / 4)) + 8 * ((M + _m) % 4 != 0);
    source = &data[tot_nodes * 8];
    target = &data[tot_nodes * 8 + M4];
    cost = &data[tot_nodes * 8 + M4 * 2];
    state = &data[tot_nodes * 8 + M4 * 3];
    memset(state, 0, M4);

    // Block size for the pivot
    _block_size = std::max(int(BLOCK_SIZE_FACTOR * std::sqrt(double(_m))),
                           MIN_BLOCK_SIZE);
    _next_arc = M0;
  }

  // Add a node balance
  void addNode(int i, int b) { balance[i] = b; }

  // Add an arc
  size_t addArc(int i, int j, int c) {
    size_t idx = M;
    source[M] = i;
    target[M] = j;
    cost[M] = c;
    state[M] = STATE_LOWER;

    // source.push_back(i);
    // target.push_back(j);
    // cost.push_back(c);
    // state.push_back(STATE_LOWER);

    M++;
    return idx;
  }

 private:
  // Initialize internal data structures
  void init() {
    // Start with the basis tree
    int tot_balance = std::accumulate(balance.begin(), balance.end(), 0);

    // Initialize artifical cost
    int ART_COST = 0;
    for (int e = M0; e < M; ++e) ART_COST = std::max(ART_COST, cost[e]);
    ART_COST = N * (ART_COST + 1);

    // Set data for the artificial root node
    root = N;
    balance[root] = -tot_balance;

    father[root] = -1;
    son[root] = 0;  // The first node is the son
    brother[root] = -1;
    pi[root] = 0;
    edge_pred[root] = -1;
    dir_pred[root] = 1;
    successors[root] = N;

    // Init spanning tree
    for (int u = 0, e = 0; u < N; ++u, ++e) {
      father[u] = root;
      son[u] = -1;
      brother[u] = u + 1;
      if (u + 1 == N) brother[u] = -1;
      successors[u] = 0;

      edge_pred[u] = e;
      state[e] = STATE_TREE;

      if (balance[u] >= 0) {
        source[e] = u;
        target[e] = root;
        cost[e] = 0;

        dir_pred[u] = DIR_UP;
        flow_pred[u] = balance[u];
        pi[u] = 0;
      } else {
        source[e] = root;
        target[e] = u;
        cost[e] = ART_COST;

        dir_pred[u] = DIR_DOWN;
        flow_pred[u] = -balance[u];
        pi[u] = ART_COST;
      }
    }
  }

  // Pricing for negative reduced cost
  bool pricing() {
    int min = 0;
    int e = _next_arc;
    int e_max = std::min(M, e + _block_size);

    while (e < M) {
      para_search(e, e_max, &state[0], &cost[0], &pi[0], &source[0], &target[0],
                  &min, &in_arc);

      if (min < 0) {
        _next_arc = in_arc;
        return true;
      }

      e = e_max;
      e_max = std::min(M, e + _block_size);
    }

    e = M0;
    e_max = std::min(_next_arc, e + _block_size);

    while (e < _next_arc) {
      para_search(e, e_max, &state[0], &cost[0], &pi[0], &source[0], &target[0],
                  &min, &in_arc);

      if (min < 0) {
        _next_arc = in_arc;
        return true;
      }

      e = e_max;
      e_max = std::min(_next_arc, e + _block_size);
    }

    return false;
  }

  bool fullPricing() {
    int best = 0;  // Tutto intero: -std::numeric_limits<float>::min();

    for (int e = M0; e < M; ++e) {
      // QUESTION: Potrei avere la lista di quelli fuori base?
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

    while (i != j) {
      if (successors[i] < successors[j]) {
        i = father[i];
      } else {
        j = father[j];
      }
    }
    join_node = i;

    // vector<bool> mark(tot_nodes, false);
    // int p = j;
    // mark[p] = true;
    // while (p != root) {
    //  p = father[p];
    //  mark[p] = true;
    //}

    // p = i;
    // while (!mark[p])
    //  p = father[p];

    // join_node = p;
  }

  // Look for leaving variable
  void findLeavingArc() {
    int i = source[in_arc];
    int j = target[in_arc];

    constexpr int MAX = std::numeric_limits<int>::max();
    delta = MAX;
    int result = 0;
    int d;

    // Search the cycle form the first node to the join node
    for (int u = i; u != join_node; u = father[u]) {
      d = flow_pred[u];
      // TODO: Questo if lo posso portare fuori ed evitare un conto e un
      // confronto
      if (dir_pred[u] == DIR_DOWN) d = MAX - d;

      if (d < delta) {
        delta = d;
        q = u;
        result = 1;
      }
    }

    // Search the cycle form the second node to the join node
    for (int u = j; u != join_node; u = father[u]) {
      d = flow_pred[u];
      // TODO: Questo if lo posso portare fuori ed evitare un conto e un
      // confronto
      if (dir_pred[u] == DIR_UP) d = MAX - d;

      if (d <= delta) {
        delta = d;
        q = u;
        result = 2;
      }
    }

    // Direction of entering arc
    if (result == 1) {
      u_in = j;
      v_in = i;
    } else {
      u_in = i;
      v_in = j;
    }
    out_arc = edge_pred[q];
    small_right = (2 * successors[q] <= N);
  }

  // Change _flow and _state vectors
  void changeFlow() {
    // para_flow(source[in_arc], join_node, &father[0], &flow_pred[0],
    //          &dir_pred[0], delta);
    // para_flow(target[in_arc], join_node, &father[0], &flow_pred[0],
    //          &dir_pred[0], -delta);

    for (int u = source[in_arc]; u != join_node; u = father[u])
      flow_pred[u] -= dir_pred[u] * delta;

    for (int u = target[in_arc]; u != join_node; u = father[u])
      flow_pred[u] += dir_pred[u] * delta;
  }

  // Update the extended API data structure
  void updateEapi() {
    int p = father[q];
    {
      // CASE 1: Update the brothers of the son of the father of
      // the leaving node
      int s = son[p];

      if (s == q) {
        son[p] = brother[q];
      } else {
        int b = brother[s];
        while (b != -1) {
          if (b == q) {
            brother[s] = brother[q];
            break;
          }
          s = b;
          b = brother[b];
        }
      }
      brother[q] = -1;
    }

    // CASE 2: Invert the path from v_in to q
    {
      if (v_in != q) {
        // First step: start with v_in (different from the other node on the
        // cycle-subpath)
        int i = v_in;
        int j = father[i];
        int h = father[j];

        // "j" becones a child of "v_in"
        father[j] = i;

        int k = son[i];
        if (k == -1)
          son[i] = j;  // O "j" diventa figlio diretto
        else {
          int l = k;  // Oppure "j" diventa un fratello del primogenito
          while (k != -1) {
            l = k;
            k = brother[k];
          }
          brother[l] = j;
        }

        // Begin the reverse the subpath
        int i0 = i;
        i = j;
        j = h;

        while (i != q) {
          h = father[j];

          father[j] = i;

          // "j" diventa uno dei figli di "i"
          k = son[i];
          if (k == i0)
            if (brother[k] == -1)
              son[i] = j;
            else {
              son[i] = brother[k];
              int l = k;
              while (k != -1) {
                l = k;
                k = brother[k];
              }
              brother[l] = j;
              brother[i0] = -1;
            }
          else {
            int l = k;
            while (k != -1) {
              l = k;
              k = brother[k];
              if (k == i0) {
                brother[l] = brother[i0];
                k = brother[k];
              }
            }
            brother[l] = j;
            brother[i0] = -1;
          }
          i0 = i;
          i = j;
          j = h;
        }

        // 3.3 Consider the neighbour of "q": remove "i0" from the sons of "q"

        // Devo aggiornare il predecessore di "q"
        if (son[q] == i0) {
          son[q] = brother[i0];
          brother[i0] = -1;
        } else {
          k = son[q];
          int l = k;
          while (k != i0) {
            l = k;
            k = brother[k];
          }
          brother[l] = brother[i0];
          brother[i0] = -1;
        }
      }

      // Update the enterinc arc
      father[v_in] = u_in;

      // "v_in" is either a direct son or a brother of the son of "u_in"
      int k = son[u_in];
      if (k == -1)
        son[u_in] = v_in;
      else {
        int l = k;
        while (k != -1) {
          l = k;
          k = brother[k];
        }
        brother[l] = v_in;
      }
    }

    {
      // 3.4: Update the successors nodes
      if (v_in != q) {
        int i = q;
        int j = father[i];

        successors[i] = successors[i] - 1 - successors[j];
        edge_pred[i] = edge_pred[j];
        flow_pred[i] = flow_pred[j];
        dir_pred[i] = -dir_pred[j];

        int i0 = i;
        i = j;
        j = father[j];

        while (i != v_in) {
          successors[i] = successors[i] + successors[i0] - successors[j];
          edge_pred[i] = edge_pred[j];
          flow_pred[i] = flow_pred[j];
          dir_pred[i] = -dir_pred[j];
          i0 = i;
          i = j;
          j = father[j];
        }

        successors[v_in] = successors[v_in] + 1 + successors[i0];
      }

      edge_pred[v_in] = in_arc;
      flow_pred[v_in] = delta;
      dir_pred[v_in] = (u_in == source[in_arc] ? DIR_DOWN : DIR_UP);

      int delta_succ = 1 + successors[v_in];

      // Update successors from "u_in" to "join_node"
      for (int u = u_in; u != join_node; u = father[u])
        successors[u] = successors[u] + delta_succ;

      // Update successors from "q" to "join_node"
      for (int u = p; u != join_node; u = father[u])
        successors[u] = successors[u] - delta_succ;

      // Update the state of the entering and leaving arcs
      state[in_arc] = STATE_TREE;
      state[out_arc] = STATE_LOWER;
    }
  }

  // Update dual variables
  void updatePotentials() {
    int dd = (pi[u_in] - dir_pred[v_in] * cost[edge_pred[v_in]]) - pi[v_in];
    if (dd == 0) return;
    if (small_right) {
      if (son[v_in] == -1) {
        pi[v_in] += dd;
        return;
      }

      int k = v_in;
      int l = 0;

      while (true) {
        do {
          do {
            pi[k] += dd;
            l = k;
            k = son[k];
          } while (k != -1);
          k = brother[l];
        } while (k != -1);

        do {
          k = father[l];
          if (k == u_in) return;
          l = k;
          k = brother[k];
        } while (k == -1);
      }
    } else {
      // The subtree containing u_in: come faccio? risalgo e poi ridiscendo?
      // updateAllPotentials();

      int k = root;
      int l = 0;

      while (true) {
        do {
          do {
            // TODO: non riesco a capire il delta??
            // pi[k] = pi[father[k]] - dir_pred[k] * cost[edge_pred[k]];
            pi[k] -= dd;
            l = k;
            k = son[k];
          } while (k != -1 && k != v_in);
          if (k == v_in) {
            l = k;
            k = brother[k];
          } else {
            k = brother[l];
            if (k == v_in) {
              l = k;
              k = brother[k];
            }
          }
        } while (k != -1);

        do {
          k = father[l];
          if (k == root) return;
          l = k;
          k = brother[k];
          if (k == v_in) {
            l = k;
            k = brother[k];
          }
        } while (k == -1);
      }
    }
  }

  // Recompute all node potentials
  void updateAllPotentials() {
    int k = root;
    pi[k] = 0;
    int l = 0;
    while (true) {
      do {
        do {
          l = k;
          k = son[k];
          if (k != -1) {
            //              fprintf(stdout, "%d -> ", k);
            pi[k] = pi[father[k]] - dir_pred[k] * cost[edge_pred[k]];
          }
        } while (k != -1);
        k = brother[l];
        if (k != -1) {
          //            fprintf(stdout, "%d -> ", k);
          pi[k] = pi[father[k]] - dir_pred[k] * cost[edge_pred[k]];
        }
      } while (k != -1);

      do {
        k = father[l];
        if (k == root) return;
        l = k;
        k = brother[k];
        if (k != -1) {
          //            fprintf(stdout, "%d -> ", k);
          pi[k] = pi[father[k]] - dir_pred[k] * cost[edge_pred[k]];
        }
      } while (k == -1);
    }
  }

  // Return the total cost of the current solution
  double totalCost() const {
    double r = 0;
    for (int u = 0; u < N; u++)
      if (source[edge_pred[u]] != root && target[edge_pred[u]] != root)
        r += cost[edge_pred[u]] * flow_pred[u];
    return r;
  }

  // return the optimal solution
  // (which format?)
  void dumpSolution() const {
    for (int u = 0; u < N; u++)
      if (flow_pred[u] > 0)
        fprintf(stdout, "{%d, %d} # %d\n", source[edge_pred[u]],
                target[edge_pred[u]], flow_pred[u]);
  }

  void dumpTable() const {
    fprintf(stdout, "\nfather :\t");
    for (int i = 0; i < tot_nodes; i++) fprintf(stdout, "%d\t", father[i]);
    fprintf(stdout, "\nson    :\t");
    for (int i = 0; i < tot_nodes; i++) fprintf(stdout, "%d\t", son[i]);
    fprintf(stdout, "\nbrother:\t");
    for (int i = 0; i < tot_nodes; i++) fprintf(stdout, "%d\t", brother[i]);
    fprintf(stdout, "\nflow_pr:\t");
    for (int i = 0; i < tot_nodes; i++) fprintf(stdout, "%d\t", flow_pred[i]);
    fprintf(stdout, "\ndir_pr :\t");
    for (int i = 0; i < tot_nodes; i++) fprintf(stdout, "%d\t", dir_pred[i]);
    fprintf(stdout, "\nedge_pr:\t");
    for (int i = 0; i < tot_nodes; i++) fprintf(stdout, "%d\t", edge_pred[i]);
    fprintf(stdout, "\nsuccess:\t");
    for (int i = 0; i < tot_nodes; i++) fprintf(stdout, "%d\t", successors[i]);
    fprintf(stdout, "\npi     :\t");
    for (int i = 0; i < tot_nodes; i++) fprintf(stdout, "%d\t", pi[i]);
    fprintf(stdout, "\n\n");
    fprintf(stdout, "Cost: %d\n", totalCost());
  }

 public:
  // Runtime in milliseconds
  double runtime() const { return _runtime; }

  double pricing_time() const { return _time_pricing; }

  // Number of iterations of simplex algorithms
  uint64_t iterations() const { return _iterations; }

  // Run the Network Simplex algorithm
  double solve() {
    auto start_tt = std::chrono::steady_clock::now();

    // Init the internal data structures
    init();

    // Loop for negative reduced cost cycles
    _time_pricing = 0.0;
    while (true) {
      _iterations++;

      // Look for negative reduced cost variables out of the basis
      auto start_t = std::chrono::steady_clock::now();
      bool stop = pricing();
      auto end_t = std::chrono::steady_clock::now();
      _time_pricing +=
          double(std::chrono::duration_cast<std::chrono::nanoseconds>(end_t -
                                                                      start_t)
                     .count()) /
          1000000000;

      if (!stop) break;

#ifdef DEBUG
      fprintf(stdout, "arc in: %d = (%d, %d)\n", in_arc, source[in_arc],
              target[in_arc]);
#endif  // DEBUG

      // Look for exiting variable
      findJoinNode();

#ifdef DEBUG
      fprintf(stdout, "join node: %d\n", join_node);
#endif  // DEBUG

      findLeavingArc();

#ifdef DEBUG
      fprintf(stdout, "arc out: %d = (%d, %d), delta: %d\n", edge_pred[q],
              source[edge_pred[q]], target[edge_pred[q]], delta);
#endif  // DEBUG

      if (delta > 0) changeFlow();

      // Update data structure
      start_t = std::chrono::steady_clock::now();
      updateEapi();
      end_t = std::chrono::steady_clock::now();
      _time_update_basis +=
          double(std::chrono::duration_cast<std::chrono::nanoseconds>(end_t -
                                                                      start_t)
                     .count()) /
          1000000000;

      start_t = std::chrono::steady_clock::now();
      updatePotentials();
      end_t = std::chrono::steady_clock::now();
      _time_update_duals +=
          double(std::chrono::duration_cast<std::chrono::nanoseconds>(end_t -
                                                                      start_t)
                     .count()) /
          1000000000;
#ifdef DEBUG
      dumpTable();
#endif  // DEBUG
    }

#ifdef DEBUG
    dumpSolution();
#endif  // DEBUG

    auto end_t = std::chrono::steady_clock::now();
    _runtime += double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_tt)
                           .count()) /
                1000;

    return totalCost();
  }
};
}  // namespace DOT
