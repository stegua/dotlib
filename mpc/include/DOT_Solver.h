/**
 * @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <numeric>

#include "DOT_EapiSimplex.h"
#include "DOT_Histogram2D.h"
#include "DOT_NetSimplex.h"
#include "DOT_NetSimplexUnit.h"

// Adde preprocessing directive
#include <ilcplex/cplex.h>

const double SCALE = 10000.0;

class graph_t {
public:
  ~graph_t(void) { free(data); }

  void reserveArcs(int _n, int _m) {
    n = _n;
    m = _m;
    idx = 0;

    size = _n + 3 * _m;
    data = (int *)malloc(size * sizeof(int));

    supply = &data[0];
    head = &data[_n + 0 * _m];
    tail = &data[_n + 1 * _m];
    cost = &data[_n + 2 * _m];
  }

  void addNode(int i, int b) { supply[i] = b; }

  void addArc(int i, int j, int c) {
    head[idx] = i;
    tail[idx] = j;
    cost[idx] = c;
    idx++;
  }

  void reset(void) { memset(data, 0, size); }

  int getSupply(int i) { return supply[i]; }

  int n;
  int m;

  int idx;
  int size;

  int *data;

  int *supply;

  int *head;
  int *tail;
  int *cost;
};

struct coprimes_t {
public:
  coprimes_t(int _v, int _w, int _c) : v(_v), w(_w), c_vw(_c) {}
  int v;
  int w;
  int c_vw;
};

typedef std::pair<int, int> int_pair;

struct pair_hash {
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2> &pair) const {
    return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
  }
};

namespace std {
template <> struct hash<std::pair<int, int>> {
  inline size_t operator()(const std::pair<int, int> &v) const {
    std::hash<int> int_hasher;
    return int_hasher(v.first) ^ int_hasher(v.second);
  }
};

} // namespace std

// Zeta block for coordinates vector
#define BLOCKSIZE 1024
#define BLOCKNUM 8

namespace DOT {

// Solver class, which wrapper the Network Simplex algorithm
class Solver {
public:
  // Standard c'tor
  Solver()
      : _runtime(0.0), _n_log(0), L(-1), verbosity(DOT_VAL_INFO), recode(""),
        opt_tolerance(0), timelimit(std::numeric_limits<double>::max()) {}

  // Setter/getter for parameters
  std::string getStrParam(const std::string &name) const {
    if (name == DOT_PAR_METHOD)
      return method;
    if (name == DOT_PAR_ALGORITHM)
      return algorithm;
    if (name == DOT_PAR_VERBOSITY)
      return verbosity;
    if (name == DOT_PAR_RECODE)
      return recode;
    return "ERROR getStrParam: wrong parameter ->" + name;
  }

  double getDblParam(const std::string &name) const {
    if (name == DOT_PAR_TIMELIMIT)
      return timelimit;
    if (name == DOT_PAR_OPTTOLERANCE)
      return opt_tolerance;
    return -1;
  }

  void setStrParam(const std::string &name, const std::string &_value) {
    std::string value(_value);
    tolower(value);

    if (name == DOT_PAR_METHOD)
      method = value;

    if (name == DOT_PAR_ALGORITHM)
      algorithm = value;

    if (name == DOT_PAR_VERBOSITY)
      verbosity = value;

    if (name == DOT_PAR_RECODE)
      recode = value;
  }

  void setDblParam(const std::string &name, double value) {
    if (name == DOT_PAR_TIMELIMIT)
      timelimit = value;

    if (name == DOT_PAR_OPTTOLERANCE)
      opt_tolerance = value;
  }

  void dumpParam() const {
    PRINT("Internal parameters: %s %s %s %s %.3f %f %s\n", method.c_str(),
          model.c_str(), algorithm.c_str(), verbosity.c_str(), timelimit,
          opt_tolerance, recode.c_str());
  }

  // Return status of the solver
  std::string status() const {
    if (_status == ProblemType::INFEASIBLE)
      return "Infeasible";
    if (_status == ProblemType::OPTIMAL)
      return "Optimal";
    if (_status == ProblemType::UNBOUNDED)
      return "Unbounded";
    if (_status == ProblemType::TIMELIMIT)
      return "TimeLimit";

    return "Undefined";
  }

  // Return runtime in milliseconds
  double runtime() const { return _runtime; }

  // Number of total iterations of simplex algorithms
  uint64_t iterations() const { return _iterations; }

  // Number of arcs in the model
  uint64_t num_arcs() const { return _num_arcs; }

  // Number of nodes in the model
  uint64_t num_nodes() const { return _num_nodes; }

  //----------------------------------------------------------------------------------------
  // Compute Kantorovich-Wasserstein distance between two measures
  std::array<double, 4> mincostflowEati(const Histogram2D &A,
                                        const Histogram2D &B,
                                        const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    // Build the graph for min cost flow
    NetSimplex simplex('F', static_cast<int>(2 * n * n),
                       static_cast<int>(n * n) * static_cast<int>(n * n));

    // Set the parameters
    simplex.setTimelimit(timelimit);
    simplex.setVerbosity(verbosity);
    simplex.setOptTolerance(0);

    auto ID = [&n](int x, int y) { return x * n + y; };

    // add first d source nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(ID(i, j), A.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(n * n + ID(i, j), -B.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
            simplex.addArc(ID(i, j), n * n + ID(i + v, j + w), p.c_vw);
          }
        }
      }

    // Init the simplex
    simplex.run();

    _runtime = simplex.runtime();
    _iterations = simplex.iterations();
    _num_arcs = simplex.num_arcs();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    double fobj = simplex.totalCost() / A.balance();

    PRINT("EAIT %s it %d UB %.6f runtime %.4f simplex %.4f "
          "num_arcs %d\n",
          msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
    ;
  }

  std::array<double, 4> mincostflowEapi(const Histogram2D &A,
                                        const Histogram2D &B,
                                        const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    // Build the graph for min cost flow
    EapiSimplex simplex(static_cast<int>(2 * n * n),
                        static_cast<int>(n * n) * static_cast<int>(n * n));

    auto ID = [&n](int x, int y) { return x * n + y; };

    // add first d source nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(ID(i, j), A.get(i, j) - B.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
            simplex.addArc(ID(i, j), ID(i + v, j + w), p.c_vw);
          }
        }
      }

    // Init the simplex
    double r = simplex.solve();
    _iterations = (int)simplex.iterations();
    _runtime = simplex.runtime();
    _num_arcs = simplex.M - simplex.M0;

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    double fobj = r / double(A.balance());

    PRINT("EAPI %s it %d UB %.6f runtime %.4f simplex %.4f "
          "num_arcs %d\n",
          msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
  }

  std::array<double, 4> bipartiteEapi(const Histogram2D &A,
                                      const Histogram2D &B,
                                      const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    double Dmax = 1;

    // Build the graph for min cost flow
    EapiSimplex simplex(static_cast<int>(2 * n * n),
                        static_cast<int>(n * n) * static_cast<int>(n * n));

    auto ID = [&n](int x, int y) { return x * n + y; };

    // add first d source nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(ID(i, j), A.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(n * n + ID(i, j), -B.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        for (int v = 0; v < n; ++v)
          for (int w = 0; w < n; ++w)
            simplex.addArc(ID(i, j), n * n + ID(v, w),
                           (i - v) * (i - v) + (w - j) * (w - j));

    // Init the simplex
    double r = simplex.solve();
    _iterations = (int)simplex.iterations();
    _runtime = simplex.runtime();
    _num_arcs = simplex.M - simplex.M0;

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    double fobj = r / double(A.balance());

    PRINT("BIP-EAPI %s it %d UB %.6f runtime %.4f simplex %.4f "
          "num_arcs %d\n",
          msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
  }

  std::array<double, 4> bipartiteEati(const Histogram2D &A,
                                      const Histogram2D &B,
                                      const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    // Build the graph for min cost flow
    NetSimplex simplex('F', static_cast<int>(2 * n * n),
                       static_cast<int>(n * n) * static_cast<int>(n * n));

    // Set the parameters
    simplex.setTimelimit(timelimit);
    simplex.setVerbosity(verbosity);
    simplex.setOptTolerance(0);

    auto ID = [&n](int x, int y) { return x * n + y; };

    // add first d source nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(ID(i, j), A.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(n * n + ID(i, j), -B.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        for (int v = 0; v < n; ++v)
          for (int w = 0; w < n; ++w)
            simplex.addArc(ID(i, j), n * n + ID(v, w),
                           (i - v) * (i - v) + (w - j) * (w - j));

    // Init the simplex
    simplex.run();

    _runtime = simplex.runtime();
    _iterations = simplex.iterations();
    _num_arcs = simplex.num_arcs();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    double fobj = simplex.totalCost() / A.balance();

    PRINT("BIP-EATI %s it %d UB %.6f runtime %.4f simplex %.4f "
          "num_arcs %d\n",
          msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
    ;
  }

  //----------------------------------------------------------------------------------------
  // Compute Kantorovich-Wasserstein distance between two measures
  double cplexNetSimplex(const Histogram2D &A, const Histogram2D &B,
                         const std::string &msg = "") {
    // Init coprimes data structure
    int n = A.getN();

    auto start_t = std::chrono::steady_clock::now();

    // Start CPLEX definition
    CPXENVptr env = NULL;
    CPXNETptr net = NULL;
    int status = 0;

    env = CPXopenCPLEX(&status);

    if (env == NULL) {
      char errmsg[CPXMESSAGEBUFSIZE];
      fprintf(stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring(env, status, errmsg);
      fprintf(stderr, "%s", errmsg);
      return -1;
    }

    CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_OFF);
    CPXsetintparam(env, CPXPARAM_Network_Tolerances_Feasibility, 1e-09);
    CPXsetintparam(env, CPXPARAM_Network_Tolerances_Optimality, 1e-09);

    net = CPXNETcreateprob(env, &status, "netCplex");

    if (net == NULL) {
      fprintf(stderr, "Failed to create network object.\n");
      fflush(stderr);
      exit(-1);
    }
    // Build Network instance
    CPXNETchgobjsen(env, net, CPX_MIN);

    vector<double> supply(2 * n * n, 0);

    auto ID = [&n](int x, int y) { return x * n + y; };

    // Add all nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        supply[ID(i, j)] = A.get(i, j);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        supply[n * n + ID(i, j)] = -B.get(i, j);

    CPXNETaddnodes(env, net, (int)supply.size(), &supply[0], NULL);

    // Add all arcs
    vector<int> tail;
    vector<int> head;
    vector<double> obj;

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
            tail.push_back(ID(i, j));
            head.push_back(n * n + ID(i + v, j + w));
            obj.push_back(p.c_vw);
          }
        }
      }
    int m = (int)tail.size();
    CPXNETaddarcs(env, net, m, &tail[0], &head[0], NULL, NULL, &obj[0], NULL);

    // Solve problem
    auto start = std::chrono::steady_clock::now();
    status = CPXNETprimopt(env, net);
    auto end = std::chrono::steady_clock::now();
    auto _runtime =
        double(
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count()) /
        1000;

    if (status) {
      fprintf(stderr, "Failed to optimize network.\n");
      fflush(stderr);
      exit(-1);
    }
    int solstat = CPXNETgetstat(env, net);

    double objval = 0.0;
    CPXNETgetobjval(env, net, &objval);

    objval = objval / double(A.balance());

    _iterations = CPXNETgetitcnt(env, net);

    /* Free up the problem as allocated by CPXNETcreateprob, if necessary */
    if (net != NULL) {
      status = CPXNETfreeprob(env, &net);
      if (status) {
        fprintf(stderr, "CPXNETfreeprob failed, error code %d.\n", status);
      }
    }

    /* Free up the CPLEX environment, if necessary */
    if (env != NULL) {
      status = CPXcloseCPLEX(&env);

      if (status) {
        char errmsg[CPXMESSAGEBUFSIZE];
        fprintf(stderr, "Could not close CPLEX environment.\n");
        CPXgeterrorstring(env, status, errmsg);
        fprintf(stderr, "%s", errmsg);
      }
    }

    // Merge with previous code
    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;
    _num_arcs = m;

    PRINT("NETCPLEX %s it %d UB %.6f runtime %.4f simplex %.4f "
          "num_arcs %d\n",
          msg.c_str(), (int)_iterations, objval, _all, _runtime,
          (int)_num_arcs);

    return _all;
  }

  double bipartiteCplex(const Histogram2D &A, const Histogram2D &B,
                        const std::string &msg = "") {
    // Init coprimes data structure
    int n = A.getN();

    auto start_t = std::chrono::steady_clock::now();

    // Start CPLEX definition
    CPXENVptr env = NULL;
    CPXNETptr net = NULL;
    int status = 0;

    env = CPXopenCPLEX(&status);

    if (env == NULL) {
      char errmsg[CPXMESSAGEBUFSIZE];
      fprintf(stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring(env, status, errmsg);
      fprintf(stderr, "%s", errmsg);
      return -1;
    }

    CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_OFF);
    CPXsetintparam(env, CPXPARAM_Network_Tolerances_Feasibility, 1e-09);
    CPXsetintparam(env, CPXPARAM_Network_Tolerances_Optimality, 1e-09);

    net = CPXNETcreateprob(env, &status, "netCplex");

    if (net == NULL) {
      fprintf(stderr, "Failed to create network object.\n");
      fflush(stderr);
      exit(-1);
    }
    // Build Network instance
    CPXNETchgobjsen(env, net, CPX_MIN);

    vector<double> supply(2 * n * n, 0);

    auto ID = [&n](int x, int y) { return x * n + y; };

    // Add all nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        supply[ID(i, j)] = A.get(i, j);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        supply[n * n + ID(i, j)] = -B.get(i, j);

    CPXNETaddnodes(env, net, (int)supply.size(), &supply[0], NULL);

    // Add all arcs
    vector<int> tail;
    vector<int> head;
    vector<double> obj;

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        for (int v = 0; v < n; ++v)
          for (int w = 0; w < n; ++w) {
            tail.push_back(ID(i, j));
            head.push_back(n * n + ID(v, w));
            obj.push_back((i - v) * (i - v) + (w - j) * (w - j));
          }

    int m = (int)tail.size();
    CPXNETaddarcs(env, net, m, &tail[0], &head[0], NULL, NULL, &obj[0], NULL);

    // Solve problem
    auto start = std::chrono::steady_clock::now();
    status = CPXNETprimopt(env, net);
    auto end = std::chrono::steady_clock::now();
    auto _runtime =
        double(
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                .count()) /
        1000;

    if (status) {
      fprintf(stderr, "Failed to optimize network.\n");
      fflush(stderr);
      exit(-1);
    }
    int solstat = CPXNETgetstat(env, net);

    double objval = 0.0;
    CPXNETgetobjval(env, net, &objval);

    objval = objval / double(A.balance());

    _iterations = CPXNETgetitcnt(env, net);

    /* Free up the problem as allocated by CPXNETcreateprob, if necessary */
    if (net != NULL) {
      status = CPXNETfreeprob(env, &net);
      if (status) {
        fprintf(stderr, "CPXNETfreeprob failed, error code %d.\n", status);
      }
    }

    /* Free up the CPLEX environment, if necessary */
    if (env != NULL) {
      status = CPXcloseCPLEX(&env);

      if (status) {
        char errmsg[CPXMESSAGEBUFSIZE];
        fprintf(stderr, "Could not close CPLEX environment.\n");
        CPXgeterrorstring(env, status, errmsg);
        fprintf(stderr, "%s", errmsg);
      }
    }

    // Merge with previous code
    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;
    _num_arcs = m;

    PRINT("BIP-PLEX %s it %d UB %.6f runtime %.4f simplex %.4f "
          "num_arcs %d\n",
          msg.c_str(), (int)_iterations, objval, _all, _runtime,
          (int)_num_arcs);

    return _all;
  }

  double colgen(const Histogram2D &A, const Histogram2D &B,
                const std::string &msg) {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    auto ID = [&n](int x, int y) { return x * n + y; };

    int N = 2 * n * n;
    vector<int> pi(N, 0);

    Vars vars(n * n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        vars[ID(i, j)].a = ID(i, j);

    Vars vnew;
    vnew.reserve(n * n);

    // Build the graph for min cost flow
    NetSimplex simplex('E', N, 0);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(ID(i, j), A.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(n * n + ID(i, j), -B.get(i, j));

    // Set the parameters
    simplex.setTimelimit(timelimit);
    simplex.setVerbosity(verbosity);
    simplex.setOptTolerance(opt_tolerance);

    _status = simplex.run();

    while (_status != ProblemType::TIMELIMIT) {
      // Take the dual values
      for (int j = 0; j < N; ++j)
        pi[j] = -simplex.potential(j);

        // Solve separation problem:
#pragma omp parallel for collapse(2)
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
          int best_v = 0;
          int best_c = -1;
          int best_n = 0; // best second node
          int h = ID(i, j);

          for (int v = 0; v < n; ++v)
            for (int w = 0; w < n; ++w) {
              int c_vw = (i - v) * (i - v) + (w - j) * (w - j);
              int violation = c_vw - pi[h] + pi[n * n + ID(v, w)];
              if (violation < best_v) {
                best_v = violation;
                best_c = c_vw;
                best_n = n * n + ID(v, w);
              }
            }

          // Store most violated cuts for element i
          vars[h].b = best_n;
          vars[h].c = best_c;
        }

      // Take all negative reduced cost variables
      vnew.clear();
      for (auto &v : vars) {
        if (v.c > -1)
          vnew.push_back(v);
        v.c = -1;
      }

      if (vnew.empty())
        break;

      std::sort(vnew.begin(), vnew.end(),
                [](const Var &v, const Var &w) { return v.c > w.c; });

      // Replace old constraints with new ones
      int new_arcs = simplex.updateArcs(vnew);

      _status = simplex.reRun();
    }

    _runtime = simplex.runtime();
    _iterations = simplex.iterations();
    _num_arcs = simplex.num_arcs();
    _num_nodes = simplex.num_nodes();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    double fobj = double(simplex.totalCost()) / A.balance();

    // Upper bound on the missed mass
    PRINT("MYCOLGEN %s it %d LB %.6f runtime %.4f simplex %.4f "
          "num_arcs %d\n",
          msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return _all;
  }

  //--------------------------------------------------------------------------
  void init_coprimes(int L, int N) {
    // La distanza la normalizzo in modo che la distanza massima sia paria  1.0
    // e poi la riscalo all'intero più vicino dopo aver moltiplicato per un
    // fattore di scala

    double Dmax = 1;

    coprimes.clear();
    for (int v = -L; v <= L; ++v)
      for (int w = -L; w <= L; ++w)
        if (!(v == 0 && w == 0) && std::gcd(v, w) == 1)
          coprimes.emplace_back(
              v, w, int(round(SCALE * sqrt(double(v * v + w * w)))));
    // Errore! non il costo quadratico, non funzinoa cosi
    coprimes.shrink_to_fit();
  }

  //--------------------------------------------------------------------------
  void init_dist_from_to(const vector<int> &data, int LO, int UP,
                         bool order = false) {
    coprimes.clear();
    for (int h = LO; h < UP; h++) {
      int L = data[h];
      for (int v = -L; v <= L; ++v)
        for (int w = -L; w <= L; ++w)
          if (pow(v, 2) + pow(w, 2) == L)
            coprimes.emplace_back(v, w, L);
    }
    if (order)
      std::sort(coprimes.begin(), coprimes.end(),
                [](const auto &a, const auto &b) { return a.c_vw < b.c_vw; });

    coprimes.shrink_to_fit();
  }

  //--------------------------------------------------------------------------
  void init_all(int n, bool order = false) {
    coprimes.clear();
    for (int v = -n + 1; v < n; ++v)
      for (int w = -n + 1; w < n; ++w)
        coprimes.emplace_back(v, w, static_cast<int>(pow(v, 2) + pow(w, 2)));
    coprimes.shrink_to_fit();

    if (order)
      std::sort(coprimes.begin(), coprimes.end(),
                [](const auto &a, const auto &b) { return a.c_vw < b.c_vw; });
  }

  //--------------------------------------------------------------------------
  void parseDimacs(const std::string &filename) {
    std::ifstream in_file(filename);

    if (!in_file) {
      fprintf(stdout, "FATAL ERROR: Cannot open file %s.\n", filename.c_str());
      exit(EXIT_FAILURE);
    }

    G.reset();

    // Read first row, and return row length
    int idx = 0;
    int m = 0;
    while (idx < m) {
      std::string line;
      std::getline(in_file, line);
      std::stringstream lineStream(line);
      std::string cell;
      // First char
      std::getline(lineStream, cell, ' ');
      if (cell == "c" || cell == "n")
        continue;
      if (cell == "p") {
        // Read "asn"
        std::getline(lineStream, cell, ' ');
        // Read "number of nodes"
        std::getline(lineStream, cell, ' ');
        int n = atoi(cell.c_str());
        // Read "number of nodes"
        std::getline(lineStream, cell, ' ');
        m = atoi(cell.c_str());

        fprintf(stdout, "read dimacs network -> n: %d, m:%d\n", n, m);
        G.reserveArcs(n, m);

        for (int i = 0; i < n / 2; ++i)
          G.addNode(i, 1);
        for (int i = n / 2; i < n; ++i)
          G.addNode(i, -1);
      }

      if (cell == "a") {
        std::getline(lineStream, cell, ' ');
        int i = atoi(cell.c_str());
        std::getline(lineStream, cell, ' ');
        int j = atoi(cell.c_str());
        std::getline(lineStream, cell, ' ');
        int c = atoi(cell.c_str());

        G.addArc(i - 1, j - 1, c);
        idx++;
      }
    }
  }

  //--------------------------------------------------------------------------
  double solveDimacs(const std::string &filename, const std::string &msg = "") {
    if (filename == "netcplex") {
      auto start_t = std::chrono::steady_clock::now();

      // Start CPLEX definition
      CPXENVptr env = NULL;
      CPXNETptr net = NULL;
      int status = 0;

      env = CPXopenCPLEX(&status);

      if (env == NULL) {
        char errmsg[CPXMESSAGEBUFSIZE];
        fprintf(stderr, "Could not open CPLEX environment.\n");
        CPXgeterrorstring(env, status, errmsg);
        fprintf(stderr, "%s", errmsg);
        return -1;
      }

      CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_OFF);
      CPXsetintparam(env, CPXPARAM_Network_Tolerances_Feasibility, 1e-09);
      CPXsetintparam(env, CPXPARAM_Network_Tolerances_Optimality, 1e-09);

      net = CPXNETcreateprob(env, &status, "netCplex");

      if (net == NULL) {
        fprintf(stderr, "Failed to create network object.\n");
        fflush(stderr);
        exit(-1);
      }
      // Build Network instance
      CPXNETchgobjsen(env, net, CPX_MIN);

      // Add all nodes
      vector<double> ss(G.n, 0);
      for (int i = 0; i < G.n; ++i)
        ss[i] = G.getSupply(i);

      CPXNETaddnodes(env, net, (int)G.n, &ss[0], NULL);

      // Add all arcs
      vector<int> tail;
      vector<int> head;
      vector<double> obj(G.m, 0);
      for (int e = 0; e < G.m; ++e)
        obj[e] = G.cost[e];

      CPXNETaddarcs(env, net, G.m, G.tail, G.head, NULL, NULL, &obj[0], NULL);

      // Solve problem
      auto start = std::chrono::steady_clock::now();
      status = CPXNETprimopt(env, net);
      auto end = std::chrono::steady_clock::now();
      auto _runtime =
          double(
              std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                  .count()) /
          1000;

      if (status) {
        fprintf(stderr, "Failed to optimize network.\n");
        fflush(stderr);
        exit(-1);
      }
      int solstat = CPXNETgetstat(env, net);

      double objval = 0.0;
      CPXNETgetobjval(env, net, &objval);

      objval = objval / double(G.n);

      _iterations = CPXNETgetitcnt(env, net);

      /* Free up the problem as allocated by CPXNETcreateprob, if necessary */
      if (net != NULL) {
        status = CPXNETfreeprob(env, &net);
        if (status) {
          fprintf(stderr, "CPXNETfreeprob failed, error code %d.\n", status);
        }
      }

      /* Free up the CPLEX environment, if necessary */
      if (env != NULL) {
        status = CPXcloseCPLEX(&env);

        if (status) {
          char errmsg[CPXMESSAGEBUFSIZE];
          fprintf(stderr, "Could not close CPLEX environment.\n");
          CPXgeterrorstring(env, status, errmsg);
          fprintf(stderr, "%s", errmsg);
        }
      }

      // Merge with previous code
      auto end_t = std::chrono::steady_clock::now();
      auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end_t - start_t)
                             .count()) /
                  1000;
      _num_arcs = G.m;

      PRINT("BIP-PLEX %s it %d UB %.6f runtime %.4f simplex %.4f "
            "num_arcs %d\n",
            msg.c_str(), (int)_iterations, objval, _all, _runtime,
            (int)_num_arcs);

      return _all;
    }
  }

  // List of pair of coprimes number between (-L, L)
  std::vector<coprimes_t> coprimes;

private:
  // Internal graph
  graph_t G;

  // Status of the solver
  ProblemType _status;

  // Runtime in milliseconds
  double _runtime;

  // Number of iterations
  uint64_t _iterations;
  uint64_t _num_nodes;
  uint64_t _num_arcs;

  // Interval for logging iterations in the simplex algorithm
  // (if _n_log=0 no logs at all)
  int _n_log;

  // Approximation parameter
  int L;

  // Method to solve the problem
  std::string method;
  // Model to solve the problem
  std::string model;
  // Algorithm to solve the corresponding problem
  std::string algorithm;
  // Verbosity of the log
  std::string verbosity;
  // Recode the coordinates as consecutive integers
  std::string recode;
  // Tolerance for pricing
  double opt_tolerance;
  // Time limit for runtime of the algorithm
  double timelimit;

}; // namespace DOT

} // namespace DOT

//		// Compute Kantorovich-Wasserstein distance between two measures
//		double tripartite(const Histogram2D& A, const Histogram2D& B) {
//			int n = A.getN();
//
//			auto start_t = std::chrono::steady_clock::now();
//
//			// Build the graph for min cost flow
//			NetSimplex simplex('F', 3 * n * n, n * n * n);
//
//			// Set the parameters
//			simplex.setTimelimit(timelimit);
//			simplex.setVerbosity(verbosity);
//			simplex.setOptTolerance(opt_tolerance);
//
//			auto ID = [&n](int x, int y) { return x * n + y; };
//
//			// add first d source nodes
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					simplex.addNode(ID(i, j), A.get(i, j));
//
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					simplex.addNode(n * n + ID(i, j), 0);
//
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					simplex.addNode(2 * n * n + ID(i, j),
//-B.get(i, j));
//
//			// First layer
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j) {
//					for (int h = 0; h < n; ++h) {
//						// fprintf(stdout, "(%d, %d)#%d
//->
//(%d, %d)#%d\t", i, j, ID(i, j), h,
//						// j,
//						//        n * n + ID(h, j));
//						simplex.addArc(ID(i, j), n * n +
// ID(h, j), (int)pow(h - i, 2));
//					}
//					//        fprintf(stdout, "\n");
//				}
//
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j) {
//					for (int h = 0; h < n; ++h) {
//						// fprintf(stdout, "(%d, %d)#%d
//->
//(%d, %d)#%d\t", i, j,
//						//        n * n + ID(i, j), i,
// h,
// 2
//*
// n
//*
// n
//+ ID(i,
// h)); 						simplex.addArc(n * n +
// ID(i, j),
// 2
// *
// n
// *
// n
// +
// ID(i, h), 							(int)pow(h - j,
// 2));
//					}
//					//       fprintf(stdout, "\n");
//				}
//
//			// Init the simplex
//			simplex.run();
//
//			_iterations = simplex.iterations();
//			_runtime = simplex.runtime();
//			_iterations = simplex.iterations();
//			_num_arcs = simplex.num_arcs();
//			_num_nodes = simplex.num_nodes();
//
//			auto end_t = std::chrono::steady_clock::now();
//			auto _all =
// double(std::chrono::duration_cast<std::chrono::milliseconds>(
// end_t - start_t) 				.count()) /
// 1000;
//
//			double fobj = simplex.totalCost() / A.balance();
//
//			PRINT("TRIPARTIE | it: %d, fobj: %.6f, runtime: %.4f
//(simplex:
//%.4f), " 				"num_arcs: %ld\n",
//_iterations, fobj, _all, _runtime, _num_arcs);
//
//			return fobj;
//		}
//
//		double tripartiteColgen(const Histogram2D& A, const Histogram2D&
// B)
//{ 			int n = A.getN(); 			auto ID =
//[&n](int x, int y) { return x * n + y; };
//
//			int N = 3 * n * n;
//			vector<int> pi(N, 0);
//
//			Vars vars(2 * n * n);
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j) {
//					vars[ID(i, j)].a = ID(i, j);
//					vars[n * n + ID(i, j)].a = n * n + ID(i,
// j);
//				}
//
//			Vars vnew;
//			vnew.reserve(2 * n * n);
//
//			auto start_t = std::chrono::steady_clock::now();
//
//			// Build the graph for min cost flow
//			NetSimplex simplex('E', 3 * n * n, 0);
//
//			// Set the parameters
//			simplex.setTimelimit(timelimit);
//			simplex.setVerbosity(verbosity);
//			simplex.setOptTolerance(opt_tolerance);
//
//			// add first d source nodes
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					simplex.addNode(ID(i, j), A.get(i, j));
//
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					simplex.addNode(n * n + ID(i, j), 0);
//
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					simplex.addNode(2 * n * n + ID(i, j),
//-B.get(i, j));
//
//			// Init the simplex
//			auto _status = simplex.run();
//
//			while (_status != ProblemType::TIMELIMIT) {
//				// Take the dual values
//				for (int j = 0; j < N; ++j)
//					pi[j] = -simplex.potential(j);
//
//				// Solve separation problem:
//#pragma omp parallel for collapse(2)
//				for (int i = 0; i < n; ++i)
//					for (int j = 0; j < n; ++j) {
//						int best_v1 = 0;
//						int best_c1 = -1;
//						int best_n1 = 0; // best second
// node 						int best_v2 = 0;
// int best_c2 = -1; int best_n2 = 0; // best second node int H1 = ID(i, j);
// int H2 = n * n
// +
// ID(i, j);
//
//						for (int h = 0; h < n; ++h) {
//							int dist = (h - i) * (h
//-
// i); 							int violation =
// dist
// - pi[H1]
// + pi[n
// *
// n
// + ID(h, j)];
//
//							if (violation < best_v1)
//{ 								best_v1 =
// violation;
// best_c1 = dist; best_n1
// =
// n
// *
// n
// + ID(h, j);
//							}
//
//							dist = (h - j) * (h -
// j); 							violation = dist
// - pi[H2]
// + pi[2
// *
// n
// *
// n
// + ID(i, h)];
//
//							if (violation < best_v2)
//{ 								best_v2 =
// violation;
// best_c2 = dist; best_n2
// =
// 2
// *
// n
// *
// n
// + ID(i, h);
//							}
//						}
//
//						vars[H1].b = best_n1;
//						vars[H1].c = best_c1;
//						vars[H2].b = best_n2;
//						vars[H2].c = best_c2;
//					}
//
//				// Take all negative reduced cost variables
//				vnew.clear();
//				for (auto& v : vars) {
//					if (v.c > -1)
//						vnew.push_back(v);
//					v.c = -1;
//				}
//
//				if (vnew.empty())
//					break;
//
//				std::sort(vnew.begin(), vnew.end(),
//					[](const Var& v, const Var& w) { return
// v.c > w.c;
//});
//
//				// Replace old constraints with new ones
//				int new_arcs = simplex.updateArcs(vnew);
//
//				_status = simplex.reRun();
//			}
//
//			_iterations = simplex.iterations();
//			_runtime = simplex.runtime();
//			_iterations = simplex.iterations();
//			_num_arcs = simplex.num_arcs();
//			_num_nodes = simplex.num_nodes();
//
//			auto end_t = std::chrono::steady_clock::now();
//			auto _all =
// double(std::chrono::duration_cast<std::chrono::milliseconds>(
// end_t - start_t) 				.count()) /
// 1000;
//
//			double fobj = simplex.totalCost() / A.balance();
//
//			PRINT("TRIPARTIE | it: %d, fobj: %.6f, runtime: %.4f
//(simplex:
//%.4f), " 				"num_arcs: %ld\n",
//_iterations, fobj, _all, _runtime, _num_arcs);
//
//			return fobj;
//		}
//
//		// Compute Kantorovich-Wasserstein distance between two measures
//		double phaseTwo(const Histogram2D& A, const Histogram2D& B) {
//			int n = A.getN();
//
//			// Compute distances
//			std::set<int> tauset;
//			for (int v = 0; v < n; ++v)
//				for (int w = 0; w < n; ++w)
//					tauset.insert(static_cast<int>(pow(v, 2)
//+ pow(w, 2)));
//
//			vector<int> tau;
//			for (auto v : tauset)
//				tau.push_back(v);
//			tau.push_back(tau.back() * 2);
//
//			fprintf(stdout, "distances: %d %d\n", tau.size(),
// tau[0]);
//
//			int idxL = 0;
//			init_dist_from_to(tau, 0, idxL);
//
//			// Build the graph for min cost flow
//			NetSimplexUnit simplex('E', static_cast<int>(2 * n * n)
//+ 1, 0);
//
//			auto ID = [&n](int x, int y) { return x * n + y; };
//			auto start_t = std::chrono::steady_clock::now();
//
//			{
//				// Set the parameters
//				simplex.setTimelimit(timelimit);
//				simplex.setVerbosity(verbosity);
//				simplex.setOptTolerance(opt_tolerance);
//
//				// add first d source nodes
//				for (int i = 0; i < n; ++i)
//					for (int j = 0; j < n; ++j)
//						simplex.addNode(ID(i, j),
// A.get(i, j));
//
//				for (int i = 0; i < n; ++i)
//					for (int j = 0; j < n; ++j)
//						simplex.addNode(n * n + ID(i,
// j), -B.get(i, j));
//
//				for (int i = 0; i < n; ++i)
//					for (int j = 0; j < n; ++j) {
//						for (const auto& p : coprimes) {
//							int v = p.v;
//							int w = p.w;
//							if (i + v >= 0 && i + v
//<
// n
//&&
// j
//+
// w
//>=
// 0
//&&
// j
//+
// w < n) {
// simplex.addArc(ID(i, j),
// n
// *
// n
// + ID(i
// + v, j + w), p.c_vw);
//							}
//						}
//					}
//
//				int it = 0;
//				int64_t fobj = 0;
//
//				// Init the simplex
//				simplex.run();
//				_iterations = simplex.iterations();
//
//				// Start separation
//				while (true) {
//					_status = simplex.reRun();
//
//					if (_status == ProblemType::TIMELIMIT)
//						break;
//
//					// Check feasibility: if feasible stop
//					auto dummyFlow = simplex.dummyFlow();
//					fprintf(stdout, "Flow: %ld, idx: %d,
// tau: %d\n", dummyFlow, idxL, tau[idxL]);
//
//					if (abs(dummyFlow) < 1e-06 || idxL >=
// tau.size()) 						break;
//
//					// Add arcs
//					init_coprimes(tau[idxL]);
//					idxL++;
//
//					for (int i = 0; i < n; ++i)
//						for (int j = 0; j < n; ++j) {
//							for (const auto& p :
// coprimes)
//{ 								int v = p.v;
// int w = p.w; if (i
// +
// v
// >=
// 0
// &&
// i
// +
// v
// <
// n
// && j + w >= 0 && j +
// w < n) {
// simplex.addArc(ID(i, j),
// n
// *
// n
// + ID(i
// + v, j + w), p.c_vw);
//								}
//							}
//						}
//
//					++it;
//				}
//
//				_runtime = simplex.runtime();
//				_iterations = simplex.iterations();
//				_num_arcs = simplex.num_arcs();
//				_num_nodes = simplex.num_nodes();
//
//				auto end_t = std::chrono::steady_clock::now();
//				auto _all =
// double(std::chrono::duration_cast<std::chrono::milliseconds>(
// end_t - start_t) 					.count()) /
// 1000;
//
//				fobj = simplex.totalCost();
//
//				PRINT("it: %d, fobj: %f, all: %f, simplex: %f,
// num_arcs:
//%ld\n", it, fobj, 					_all, _runtime,
//_num_arcs);
//
//				return fobj;
//			}
//
//			// ------------------------ PHASE TWO
//			// -----------------------------------
//			///*idxL = 3;
//			// init_coprimes(tau[idxL]);
//			// idxL++;*/
//
//			// idxL = 1;
//			// init_dist_upto(tau[idxL]);
//			// idxL++;
//
//			// fprintf(stdout, "distances: %d %d\n", tau.size(),
// tau[idxL
//- 1]);
//
//			//// Build the graph for min cost flow
//			// NetSimplex<int, int> simplexTwo('E',
// static_cast<int>(2
//*
// n
//* n + 1),
//			//	0);
//
//			// for (size_t i = 0; i < n; ++i)
//			//	for (size_t j = 0; j < n; ++j)
//			//		simplexTwo.addNode(ID(i, j), A.get(i,
// j));
//
//			// for (size_t i = 0; i < n; ++i)
//			//	for (size_t j = 0; j < n; ++j)
//			//		simplexTwo.addNode(n * n + ID(i, j),
//-B.get(i, j));
//
//			// simplexTwo.addNode(2 * n * n, 0);
//
//			// for (size_t i = 0; i < n; ++i)
//			//	for (size_t j = 0; j < n; ++j) {
//			//		for (const auto& p : coprimes) {
//			//			int v = p.v;
//			//			int w = p.w;
//			//			if (i + v >= 0 && i + v < n && j
//+
// w
//>=
// 0
//&&
// j
//+
//			// w
//			//<
//			// n) {
// simplexTwo.addArc(ID(i, j),
// n
//*
//			// n
//			// + ID(i
//			//+ v,
//			// j
//			//+ w), p.c_vw);
//			//			}
//			//		}
//			//	}
//
//			NetSimplex simplexTwo(simplex);
//
//			// Arcs to be updated
//			vector<size_t> dummy_arcs;
//			dummy_arcs.reserve(n * n);
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					dummy_arcs.push_back(simplexTwo.addArc(ID(i,
// j),
// 2
//*
// n
//* n, tau[idxL]));
//
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					simplexTwo.addArc(2 * n * n, ID(i, j),
// 0);
//
//			simplexTwo.recomputePotential();
//
//			// Set the parameters
//			simplexTwo.setTimelimit(timelimit);
//			simplexTwo.setVerbosity(verbosity);
//			simplexTwo.setOptTolerance(opt_tolerance);
//
//			_status = simplexTwo.reRun();
//
//			double totF = A.balance();
//			while (_status != ProblemType::TIMELIMIT) {
//				// Check feasibility: if feasible stop
//				auto dummyFlow =
// simplexTwo.computeDummyFlow(dummy_arcs);
// fprintf(stdout, "Flow: %d, idx: %d, tau: %d, fobj: %.5f\n", dummyFlow,
// idxL, tau[idxL - 1], double(dummyFlow) / totF);
//
//				if (dummyFlow < 1 || idxL >= tau.size() - 1)
//					break;
//
//				// Add arcs
//				init_coprimes(tau[idxL]);
//				idxL++;
//
//				for (int i = 0; i < n; ++i)
//					for (int j = 0; j < n; ++j) {
//						for (const auto& p : coprimes) {
//							int v = p.v;
//							int w = p.w;
//							if (i + v >= 0 && i + v
//<
// n
//&&
// j
//+
// w
//>=
// 0
//&&
// j
//+
// w < n) {
// simplexTwo.addArc(ID(i, j),
// n
// *
// n
// + ID(i
// + v, j + w), p.c_vw);
//							}
//						}
//					}
//
//				simplexTwo.updateArcs(dummy_arcs, tau[idxL]);
//				simplexTwo.recomputePotential();
//
//				_status = simplexTwo.reRun();
//			}
//
//			_runtime = simplexTwo.runtime();
//			_iterations += simplexTwo.iterations();
//			_num_arcs = simplexTwo.num_arcs();
//			_num_nodes = simplexTwo.num_nodes();
//
//			auto end_t = std::chrono::steady_clock::now();
//			auto _all =
// double(std::chrono::duration_cast<std::chrono::milliseconds>(
// end_t - start_t) 				.count()) /
// 1000;
//
//			double fobj = double(simplexTwo.totalCost()) /
// A.balance();
//
//			PRINT("NEARBY    | it: %d, fobj: %.6f, runtime: %.4f
//(simplex:
//%.4f), " 				"num_arcs: %ld\n",
//_iterations, fobj, _all, _runtime, _num_arcs);
//
//			return fobj;
//		}
//
//		// Compute Kantorovich-Wasserstein distance between two measures
//		double phaseOne(const Histogram2D& A, const Histogram2D& B) {
//			int n = A.getN();
//
//			// Compute distances
//			std::set<int> tauset;
//			for (int v = 0; v < n; ++v)
//				for (int w = 0; w < n; ++w)
//					tauset.insert(static_cast<int>(pow(v, 2)
//+ pow(w, 2)));
//
//			vector<int> tau;
//			for (auto v : tauset)
//				tau.push_back(v);
//
//			tau.push_back(tau.back() * 2);
//
//			fprintf(stdout, "distances: %ld %dl\n", tau.size(),
// tau[0]);
//
//			auto ID = [&n](int x, int y) { return x * n + y; };
//
//			int N = 2 * n * n;
//			vector<int> pi(N, 0);
//
//			Vars vars(N);
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					vars[ID(i, j)].a = ID(i, j);
//
//			Vars vnew;
//			vnew.reserve(N);
//
//			int TT = 5;
//			int idxLO = 0;
//			int idxUP = tau.size() / 10;
//			init_dist_from_to(tau, idxLO, idxUP);
//
//			// Build the graph for min cost flow
//			NetSimplexUnit simplex('E', static_cast<int>(2 * n * n),
// 0);
//
//			auto start_t = std::chrono::steady_clock::now();
//
//			// Set the parameters
//			simplex.setTimelimit(timelimit);
//			simplex.setVerbosity(verbosity);
//			simplex.setOptTolerance(opt_tolerance);
//
//			// add first d source nodes
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					simplex.addNode(ID(i, j), A.get(i, j));
//
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					simplex.addNode(n * n + ID(i, j),
//-B.get(i, j));
//
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j) {
//					for (const auto& p : coprimes) {
//						int v = p.v;
//						int w = p.w;
//						if (i + v >= 0 && i + v < n && j
//+
// w
//>=
// 0
//&&
// j
//+
// w
//<
// n) { simplex.addArc(ID(i, j),
// n
// *
// n
// + ID(i
// + v,
// j
// + w), p.c_vw);
//						}
//					}
//				}
//
//			int64_t fobj = 0;
//
//			// Init the simplex
//			simplex.run();
//
//			while (false && _status != ProblemType::TIMELIMIT) {
//				// Check feasibility: if feasible stop
//				auto dummyFlow = simplex.dummyFlow();
//				fprintf(stdout, "Flow: %ld, idx: %d, tau:
//[%d,%d)\n", dummyFlow, idxUP, tau[idxLO], tau[idxUP]);
//
//				if (dummyFlow == 0 || idxUP >= tau.size() - 1)
//					break;
//
//				// Add arcs
//				idxLO = idxUP;
//				idxUP = std::min(idxUP + TT, (int)tau.size());
//				init_dist_from_to(tau, idxLO, idxUP);
//
//				for (int i = 0; i < n; ++i)
//					for (int j = 0; j < n; ++j) {
//						for (const auto& p : coprimes) {
//							int v = p.v;
//							int w = p.w;
//							if (i + v >= 0 && i + v
//<
// n
//&&
// j
//+
// w
//>=
// 0
//&&
// j
//+
// w < n) {
// simplex.addArc(ID(i, j),
// n
// *
// n
// + ID(i
// + v, j + w), p.c_vw);
//							}
//						}
//					}
//				_status = simplex.reRun();
//			}
//
//			_runtime = simplex.runtime();
//			_iterations = simplex.iterations();
//			_num_arcs = simplex.num_arcs();
//			_num_nodes = simplex.num_nodes();
//
//			auto end_t = std::chrono::steady_clock::now();
//			auto _all =
// double(std::chrono::duration_cast<std::chrono::milliseconds>(
// end_t - start_t) 				.count()) /
// 1000;
//
//			fobj = simplex.totalCost();
//
//			PRINT("PHASE ONE  | it: %d, fobj: %.6f, runtime: %.4f
//(simplex:
//%.4f), " 				"num_arcs: %ld\n",
//_iterations, fobj, _all, _runtime, _num_arcs);
//
//			return fobj;
//		}
//
//		double colgenOld(const Histogram2D& A, const Histogram2D& B, int
// idxL, 			const std::string& msg) { int n = A.getN();
//
//			// Compute distances
//			std::set<int> tauset;
//			for (int v = 0; v < n; ++v)
//				for (int w = 0; w < n; ++w)
//					tauset.insert(static_cast<int>(pow(v, 2)
//+ pow(w, 2)));
//
//			vector<int> tau;
//			for (auto v : tauset)
//				tau.push_back(v);
//			fprintf(stdout, "distances: %d\n", tau.size());
//
//			//    tau.push_back(tau.back() * 2);
//
//			//    size_t idxL =;
//			int TT = std::min<int>(1024, static_cast<int>(tau.size()
//- 1)); 			init_dist_from_to(tau, 0, TT);
//
//			fprintf(stdout, "coprimes size: %d\n", coprimes.size());
//
//			//coprimes.resize(128);
//
//			auto ID = [&n](int x, int y) { return x * n + y; };
//
//			int N = 2 * n * n;
//			vector<int> pi(N + 1, 0);
//
//			Vars vars(N + 1);
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					vars[ID(i, j)].a = ID(i, j);
//			vars[N].a = N;
//
//			Vars vnew;
//			vnew.reserve(N + 1);
//
//			auto start_t = std::chrono::steady_clock::now();
//
//			// Build the graph for min cost flow
//			NetSimplex simplex('E', N + 1, 0);
//
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					simplex.addNode(ID(i, j), A.get(i, j));
//
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					simplex.addNode(n * n + ID(i, j),
//-B.get(i, j));
//
//			simplex.addNode(N, 0);
//
//			vector<size_t> left_arcs, right_arcs;
//			left_arcs.reserve(n * n);
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					left_arcs.push_back(simplex.addArc(ID(i,
// j), N, tau[TT]));
//
//			for (int i = 0; i < n; ++i)
//				for (int j = 0; j < n; ++j)
//					right_arcs.push_back(simplex.addArc(N, n
//*
// n
//+ ID(i, j), 0));
//
//			// Set the parameters
//			simplex.setTimelimit(timelimit);
//			simplex.setVerbosity(verbosity);
//			simplex.setOptTolerance(opt_tolerance);
//
//			_status = simplex.run();
//
//			while (_status != ProblemType::TIMELIMIT) {
//				// Take the dual values
//				for (int j = 0; j < N; ++j)
//					pi[j] = -simplex.potential(j);
//
//				// Solve separation problem:
////#pragma omp parallel for collapse(2)
//				for (int i = 0; i < n; ++i)
//					for (int j = 0; j < n; ++j) {
//						int best_v = 0;
//						int best_c = -1;
//						int best_n = 0; // best second
// node 						int h = ID(i, j);
//						//for (const auto& p : coprimes)
//						for (int PP = 0; PP <
// std::min<int>(coprimes.size(), BLOCKSIZE * BLOCKNUM); PP++)
//						{
//							const auto& p =
// coprimes[PP]; 							int v =
// p.v; int
// w
// = p.w; 							if (i + v >= 0
// &&
// i
// +
// v
// <
// n
// &&
// j
// +
// w
//>= 0 && j + w < n) {
// int violation = p.c_vw
//- pi[h] + pi[n * n + ID(i + v, j
//+ w)]; 								if
//(violation < best_v)
//{
// best_v = violation;
// best_c = p.c_vw; best_n = n * n + ID(i + v, j + w);
//								}
//							}
//						}
//
//						// Store most violated cuts for
// element i 						vars[h].b = best_n;
// vars[h].c = best_c;
//					}
//
//				// Take all negative reduced cost variables
//				vnew.clear();
//				for (auto& v : vars) {
//					fprintf(stdout, "%d %d\t", v.b, v.c);
//					if (v.c > -1)
//						vnew.push_back(v);
//					v.c = -1;
//				}
//				fprintf(stdout, "\n");
//				fflush(stdout);
//
//				if (vnew.empty())
//					break;
//
//				std::sort(vnew.begin(), vnew.end(),
//					[](const Var& v, const Var& w) { return
// v.c > w.c;
//});
//
//				// Replace old constraints with new ones
//				int new_arcs = simplex.addArcs(vnew);
//
//				_status = simplex.reRun();
//			}
//
//			_runtime = simplex.runtime();
//			_iterations = simplex.iterations();
//			_num_arcs = simplex.num_arcs();
//			_num_nodes = simplex.num_nodes();
//
//			auto end_t = std::chrono::steady_clock::now();
//			auto _all =
// double(std::chrono::duration_cast<std::chrono::milliseconds>(
// end_t - start_t) 				.count()) /
// 1000;
//
//			double fobj = double(simplex.totalCost()) / A.balance();
//
//			// Upper bound on the missed mass
//			auto unmoved = 0.0;//
// double(simplex.computeDummyFlow(left_arcs)) / A.balance();
// double delta = 0.0;
//
//			vector<double> Aflow(n * n, 0.0);
//			vector<double> Bflow(n * n, 0.0);
//
//			Histogram2D AA(A);
//			Histogram2D BB(B);
//
//			if (unmoved > 0) {
//				for (int i = 0; i < n; ++i)
//					for (int j = 0; j < n; ++j) {
//						Aflow[ID(i, j)] =
// fabs(simplex.arcFlow(left_arcs[ID(i, j)]));
// AA.set(i, j, Aflow[ID(i, j)]); Aflow[ID(i, j)] = Aflow[ID(i, j)] /
// A.balance();
//
//						Bflow[ID(i, j)] =
// fabs(simplex.arcFlow(right_arcs[ID(i, j)]));
// BB.set(i, j, Bflow[ID(i, j)]); Bflow[ID(i, j)] = Bflow[ID(i, j)] /
// A.balance();
//					}
//
//				for (int i = 0; i < n; ++i)
//					for (int j = 0; j < n; ++j)
//						for (int v = 0; v < n; ++v)
//							for (int w = 0; w < n;
//++w) 								delta +=
// std::max(0, static_cast<int>(pow(i - v, 2) + pow(j - w, 2)) -
// tau[TT + 1])
//* 								(Aflow[ID(i,
// j)])
//* (Bflow[ID(v, w)])
/// unmoved;
//			}
//
//			// delta = findUB(AA, BB) / A.balance();
//
//			PRINT("COLGEN %s it %lld LB %.6f UB %.6f runtime %.4f
// simplex
//%.4f " 				"num_arcs %lld idx %d tau %d maxtau %d
// residual
//%.6f\n", 				msg.c_str(), _iterations, fobj, delta,
//_all, _runtime, _num_arcs, TT, 				tau[TT],
// tau[tau.size() - 1], unmoved);
//
//			return fobj;
//		}