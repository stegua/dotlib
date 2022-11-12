/**
 * @fileoverview Copyright (c) 2019-2022, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <cstring>

#include "DOT_Commons.h"

// Gurobi Simplex
extern "C" {
#include "gurobi_c.h"
}

namespace DOT {

class SolverGurobi {
 public:
  SolverGurobi(void) {}

  // Free all pointers
  ~SolverGurobi(void) {
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
    supply = &data[3 * _m];
    memset(supply, 0.0, _n * sizeof(double));
  }

  // Reset memory for reusing the allocated memory
  void reset(void) {
    memset(arcs, 0, size_arcs);
    memset(data, 0, size_data);
  }

  // Add node
  void addNode(size_t i, double b) { supply[i] = b; }

  // Add arc
  inline void addArc(size_t i, size_t j, double _cost, double _capLB,
                     double _capUB) noexcept {
    head[idx] = i;
    tail[idx] = j;

    cost[idx] = _cost;
    capLB[idx] = _capLB;
    capUB[idx] = _capUB;
    idx++;
  }

  // Solve the instance
  double solve(const std::string& msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    // Start CPLEX definition
    GRBenv* env = NULL;
    GRBmodel* model = NULL;
    int status = 0;

    GRBloadenv(&env, NULL);
    GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
    GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
    GRBsetdblparam(env, GRB_DBL_PAR_FEASIBILITYTOL, 1e-09);
    GRBsetdblparam(env, GRB_DBL_PAR_OPTIMALITYTOL, 1e-09);

    // Arrays for passing data to Gurobi
    char* le = (char*)malloc(sizeof(char) * (n + 1));
    int* ind = (int*)malloc(sizeof(int) * (n + 1));
    double* val = (double*)malloc(sizeof(double) * (n + 1));
    double* xbar = (double*)malloc(sizeof(double) * m);
    double* pert = (double*)malloc(sizeof(double) * n);

    // Variables and objective function
    GRBnewmodel(env, &model, "ot", m, cost, NULL, NULL, NULL, NULL);

    // Add constraints
    //    for (int i = 0; i < n; i++)

    // Solve problem
    auto start = std::chrono::steady_clock::now();
    status = GRBoptimize(model);
    auto end = std::chrono::steady_clock::now();
    auto _runtime =
        double(
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count()) /
        1000;

    if (status) {
      fprintf(stderr, "Failed to optimize network. Status: %d\n", status);
      fflush(stderr);
      exit(-1);
    }

    // fprintf(stdout,
    //        "BIP-GURO %s it %d UB %.6f runtime %.4f simplex %.4f "
    //        "num_arcs %d\n",
    //        msg.c_str(), (int)_iterations, objval, _all, _runtime, (int)m);

    return _all;
  }

  // Dump the graph to a file
  void dump(void) const {
    for (size_t u = 0; u < n; u++) fprintf(stdout, "%ld %f\n", u, supply[u]);

    for (size_t e = 0; e < m; e++)
      fprintf(stdout,
              "(head: %ld, tail: %ld, lb: %.3f, ub: %.3f, cost: %.3f)\n",
              (int)head[e], (int)tail[e], capLB[e], capUB[e], cost[e]);
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

void parseDIMACS(const std::string filename, size_t& len,
                 SolverGurobi& solver) noexcept {
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
      buf = std::strchr(buf, sep);
      buf = std::strchr(buf, eol);
    } else {
      if (buf[0] == 'p') {
        buf = std::strchr(buf, sep);
        ++buf;
        buf = std::strchr(buf, sep);
        ++buf;
        solver.n = std::atol(buf);
        buf = std::strchr(buf, sep);
        ++buf;
        solver.m = std::atol(buf);
        buf = std::strchr(buf, eol);
        ++buf;
        fprintf(stdout, "%d %d\n", solver.n, solver.m);
        break;
      }
    }
  }

  // check
  if (solver.n > 0 && solver.m > 0) solver.reserveArcs(solver.n, solver.m);

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

    solver.addNode(_idx - 1, -_s);
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

    solver.addArc(i - 1, j - 1, _c, _lb, _ub);

    buf = std::strchr(buf, eol);
    ++buf;
  }
}

}  // namespace DOT
