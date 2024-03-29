/**
 * @fileoverview Copyright (c) 2019-2022, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <ilcplex/cplex.h>

#include <algorithm>
#include <numeric>

#include "DOT_EapiSimplex.h"
#include "DOT_Histogram2D.h"
#include "DOT_NetSimplex.h"
#include "DOT_NetSimplexUnit.h"

#define LEMON_ONLY_TEMPLATES
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>

// Gurobi Simplex
extern "C" {
#include "gurobi_c.h"
}

#include "yocta_logger.hh"

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

  void dump(void) const {
    for (int u = 0; u < n; u++) fprintf(stdout, "%d %d\n", u, supply[u]);
    for (int e = 0; e < m; e++)
      fprintf(stdout, "(%d, %d)#%d\n", head[e], tail[e], cost[e]);
  }

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
template <>
struct hash<std::pair<int, int>> {
  inline size_t operator()(const std::pair<int, int> &v) const {
    std::hash<int> int_hasher;
    return int_hasher(v.first) ^ int_hasher(v.second);
  }
};

}  // namespace std

// Zeta block for coordinates vector
#define BLOCKSIZE 1024
#define BLOCKNUM 8

namespace DOT {
yocta::Logger logger;

// Solver class, which wrapper the Network Simplex algorithm
class Solver {
 public:
  // Standard c'tor
  Solver()
      : _runtime(0.0),
        _n_log(0),
        L(-1),
        verbosity(DOT_VAL_INFO),
        recode(""),
        opt_tolerance(0),
        timelimit(std::numeric_limits<double>::max()) {}

  // Setter/getter for parameters
  std::string getStrParam(const std::string &name) const {
    if (name == DOT_PAR_METHOD) return method;
    if (name == DOT_PAR_ALGORITHM) return algorithm;
    if (name == DOT_PAR_VERBOSITY) return verbosity;
    if (name == DOT_PAR_RECODE) return recode;
    return "ERROR getStrParam: wrong parameter ->" + name;
  }

  double getDblParam(const std::string &name) const {
    if (name == DOT_PAR_TIMELIMIT) return timelimit;
    if (name == DOT_PAR_OPTTOLERANCE) return opt_tolerance;
    return -1;
  }

  void setStrParam(const std::string &name, const std::string &_value) {
    std::string value(_value);
    tolower(value);

    if (name == DOT_PAR_METHOD) method = value;

    if (name == DOT_PAR_ALGORITHM) algorithm = value;

    if (name == DOT_PAR_VERBOSITY) verbosity = value;

    if (name == DOT_PAR_RECODE) recode = value;
  }

  void setDblParam(const std::string &name, double value) {
    if (name == DOT_PAR_TIMELIMIT) timelimit = value;

    if (name == DOT_PAR_OPTTOLERANCE) opt_tolerance = value;
  }

  void dumpParam() const {
    PRINT("Internal parameters: %s %s %s %s %.3f %f %s\n", method.c_str(),
          model.c_str(), algorithm.c_str(), verbosity.c_str(), timelimit,
          opt_tolerance, recode.c_str());
  }

  // Return status of the solver
  std::string status() const {
    if (_status == ProblemType::INFEASIBLE) return "Infeasible";
    if (_status == ProblemType::OPTIMAL) return "Optimal";
    if (_status == ProblemType::UNBOUNDED) return "Unbounded";
    if (_status == ProblemType::TIMELIMIT) return "TimeLimit";

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
      for (int j = 0; j < n; ++j) simplex.addNode(ID(i, j), A.get(i, j));

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

    PRINT(
        "EAIT %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d\n",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
    ;
  }

  //--------------------------------------------------------------------------
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

    PRINT(
        "EAPI %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d\n",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
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
    CPXsetdblparam(env, CPXPARAM_Network_Tolerances_Feasibility, 1e-09);
    CPXsetdblparam(env, CPXPARAM_Network_Tolerances_Optimality, 1e-09);

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
      for (int j = 0; j < n; ++j) supply[ID(i, j)] = A.get(i, j);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) supply[n * n + ID(i, j)] = -B.get(i, j);

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

    DOT::logger.info(
        "NETCPLEX %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, objval, _all, _runtime, (int)_num_arcs);

    return _all;
  }

  //--------------------------------------------------------------------------
  std::array<double, 4> bipartiteEapi(const Histogram2D &A,
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
      for (int j = 0; j < n; ++j) simplex.addNode(ID(i, j), A.get(i, j));

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

    DOT::logger.info(
        "BIP-EAPI %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
  }

  //--------------------------------------------------------------------------
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
      for (int j = 0; j < n; ++j) simplex.addNode(ID(i, j), A.get(i, j));

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

    DOT::logger.info(
        "BIP-EATI %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
    ;
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
    CPXsetdblparam(env, CPXPARAM_Network_Tolerances_Feasibility, 1e-09);
    CPXsetdblparam(env, CPXPARAM_Network_Tolerances_Optimality, 1e-09);

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
      for (int j = 0; j < n; ++j) supply[ID(i, j)] = A.get(i, j);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) supply[n * n + ID(i, j)] = -B.get(i, j);

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

    DOT::logger.info(
        "BIP-PLEX %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, objval, _all, _runtime, (int)_num_arcs);

    return _all;
  }

  double bipartiteGurobi(const Histogram2D &A, const Histogram2D &B,
                         const std::string &msg = "", int method = 1) {
    // Init coprimes data structure
    int n = A.getN();

    auto start_t = std::chrono::steady_clock::now();

    GRBenv *env = NULL;
    GRBmodel *model = NULL;
    int status = 0;

    GRBloadenv(&env, NULL);
    GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
    if (method == 1) GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
    if (method == 2) {
      GRBsetintparam(env, GRB_INT_PAR_CROSSOVER, 0);
      GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER);
    }
    if (method == 3) GRBsetintparam(env, GRB_INT_PAR_NETWORKALG, 1);

    // GRBsetdblparam(env, GRB_DBL_PAR_FEASIBILITYTOL, 1e-09);
    // GRBsetdblparam(env, GRB_DBL_PAR_OPTIMALITYTOL, 1e-09);

    auto ID = [&n](int x, int y) { return x * n + y; };

    vector<double> obj;

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        for (int v = 0; v < n; ++v)
          for (int w = 0; w < n; ++w) {
            obj.push_back((i - v) * (i - v) + (w - j) * (w - j));
          }

    int m = (int)obj.size();

    // Variables and objective function
    GRBnewmodel(env, &model, "ot", m, &obj[0], NULL, NULL, NULL, NULL);

    // Add constraints
    vector<double> val(n * n, 1.0);
    vector<int> ind;
    ind.reserve(n * n);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        ind.clear();
        for (int v = 0; v < n; ++v)
          for (int w = 0; w < n; ++w)
            ind.push_back(ID(i, j) * n * n + ID(v, w));

        GRBaddconstr(model, n * n, &ind[0], &val[0], GRB_EQUAL, A.get(i, j),
                     NULL);
      }

    for (int v = 0; v < n; ++v)
      for (int w = 0; w < n; ++w) {
        ind.clear();
        for (int i = 0; i < n; ++i)
          for (int j = 0; j < n; ++j)
            ind.push_back(ID(i, j) * n * n + ID(v, w));

        GRBaddconstr(model, n * n, &ind[0], &val[0], GRB_EQUAL, B.get(v, w),
                     NULL);
      }

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
      fprintf(stderr, "Failed to optimize network.\n");
      fflush(stderr);
      exit(-1);
    }

    double objval = 0.0;
    GRBgetdblattr(model, "ObjVal", &objval);

    objval = objval / double(A.balance());

    double iters = 0.0;
    GRBgetdblattr(model, "IterCount", &iters);
    _iterations = iters;

    if (model != NULL) GRBfreemodel(model);

    if (env != NULL) GRBfreeenv(env);

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;
    _num_arcs = m;

    DOT::logger.info(
        "BIP-GUR%d %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        method, msg.c_str(), (int)_iterations, objval, _all, _runtime,
        (int)_num_arcs);

    return _all;
  }

  double bipartiteColgen(const Histogram2D &A, const Histogram2D &B,
                         const std::string &msg) {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    auto ID = [&n](int x, int y) { return x * n + y; };

    int N = 2 * n * n;
    vector<int> pi(N, 0);
    pi.shrink_to_fit();

    int N2 = n * n;
    vector<int> X(N2, 0);
    vector<int> I(N2, 0);
    vector<int> J(N2, 0);
    X.shrink_to_fit();
    I.shrink_to_fit();
    J.shrink_to_fit();

    Vars vars(n * n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        auto idx = ID(i, j);
        X[idx] = i * i + j * j;
        I[idx] = i;
        J[idx] = j;
        vars[idx].a = idx;
      }

    Vars vnew;
    vnew.reserve(n * n);

    // Build the graph for min cost flow
    NetSimplex simplex('E', N, 0);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) simplex.addNode(ID(i, j), A.get(i, j));

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
      for (int j = 0; j < N; ++j) pi[j] = -simplex.potential(j);

        // Solve separation problem
#pragma omp parallel for collapse(2)
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
          int best_v = 0;
          int best_c = -1;
          int best_n = 0;  // best second node
          int h = ID(i, j);

          // i^2 + j^2 + w^2 + v^2 - 2*i*v - 2*w*j
          int tmp = i * i + j * j - pi[h];
          for (int v = 0; v < n; ++v) {
            int tmp0 = tmp + v * v - 2 * i * v;
            for (int w = 0; w < n; ++w) {
              // int c_vw = (i - v) * (i - v) + (w - j) * (w - j);
              //   ==>  i^2 + j^2 + w^2 + v^2 - 2*i*v - 2*w*j
              int c_vw = tmp0 + w * w - 2 * w * j;
              // int violation = c_vw - pi[h] + pi[n2 + ID(v, w)];
              int violation = c_vw + pi[N2 + ID(v, w)];
              if (violation < best_v) {
                best_v = violation;
                best_c = c_vw + pi[h];
                best_n = N2 + ID(v, w);
              }
            }
          }

          // Store most violated cuts for element i
          vars[h].b = best_n;
          vars[h].c = best_c;
        }

      // Take all negative reduced cost variables
      vnew.clear();
      for (auto &v : vars) {
        if (v.c > -1) vnew.push_back(v);
        v.c = -1;
      }

      if (vnew.empty()) break;

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
    DOT::logger.info(
        "MYCOLGEN %s it %d LB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return _all;
  }

  double bipartiteLemon(const Histogram2D &A, const Histogram2D &B,
                        const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    using namespace lemon;
    typedef lemon::ListDigraph Graph;
    Graph g;

    auto ID = [&n](int x, int y) { return x * n + y; };

    int d = n * n;

    // add d nodes for each histrogam (d+1) source, (d+2) target
    std::vector<Graph::Node> nodes_A;
    nodes_A.reserve(d);
    // add first d source nodes for first partition
    for (size_t i = 0; i < d; ++i) nodes_A.emplace_back(g.addNode());

    std::vector<Graph::Node> nodes_B;
    nodes_B.reserve(d);
    // add first d source nodes for second partition
    for (size_t i = 0; i < d; ++i) nodes_B.emplace_back(g.addNode());

    // Add arcs for complete bipartite graph
    std::vector<Graph::Arc> arcs;
    arcs.reserve(d * d);
    std::vector<int64_t> arcs_costs;
    arcs_costs.reserve(d * d);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        for (int v = 0; v < n; ++v)
          for (int w = 0; w < n; ++w) {
            arcs.emplace_back(g.addArc(nodes_A[ID(i, j)], nodes_B[ID(v, w)]));
            arcs_costs.emplace_back(pow((i - v), 2) + pow((j - w), 2));
          }
      }

    // lower and upper bounds, cost
    ListDigraph::ArcMap<int64_t> l_i(g), u_i(g);
    ListDigraph::ArcMap<int64_t> c_i(g);

    // FLow balance
    ListDigraph::NodeMap<int64_t> b_i(g);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        b_i[nodes_A[ID(i, j)]] = +A.get(i, j);
        b_i[nodes_B[ID(i, j)]] = -B.get(i, j);
      }

    // Add all edges
    for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto &a = arcs[i];
      l_i[a] = 0;
      u_i[a] = std::numeric_limits<int64_t>::max() - 1;
      c_i[a] = arcs_costs[i];
    }

    // Build the graph for min cost flow
    lemon::NetworkSimplex<Graph, int64_t, int64_t> simplex(g);

    // set lower/upper bounds, cost
    simplex.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

    NetworkSimplex<Graph, int64_t, int64_t>::ProblemType ret = simplex.run();

    switch (ret) {
      case NetworkSimplex<Graph>::INFEASIBLE:
        fprintf(stdout, "INFEASIBLE");
        break;
      case NetworkSimplex<Graph>::UNBOUNDED:
        fprintf(stdout, "UNBOUNDED");
        break;
    }

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    double fobj = double(simplex.totalCost()) / A.balance();

    DOT::logger.info(
        "BIP-LEMO %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return _all;
  }

  //--------------------------------------------------------------------------
  void init_coprimes(int L, int N) {
    // La distanza la normalizzo in modo che la distanza massima sia
    // paria  1.0 e poi la riscalo all'intero pi� vicino dopo aver
    // moltiplicato per un fattore di scala

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
          if (pow(v, 2) + pow(w, 2) == L) coprimes.emplace_back(v, w, L);
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
    int m = std::numeric_limits<int>::max();
    while (idx < m) {
      std::string line;
      std::getline(in_file, line);
      std::stringstream lineStream(line);
      std::string cell;
      // First char
      std::getline(lineStream, cell, ' ');

      if (cell == "c" || cell == "n") continue;
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

        for (int i = 0; i < n / 2; ++i) G.addNode(i, 1);
        for (int i = n / 2; i < n; ++i) G.addNode(i, -1);
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
    fflush(stdout);
  }

  //--------------------------------------------------------------------------
  double solveDimacs(const std::string &filename, const std::string &msg = "") {
    //#ifdef __MY_CUDA
    if (filename == "colgen_cuda") {
#ifdef __MY_CUDA
      auto start_t = std::chrono::steady_clock::now();

      int N = G.n;

      // Dual multipliers
      size_t Nm = G.m * sizeof(int);
      int *d_cost;
      int *d_head;
      int *d_tail;
      cudaMalloc((void **)&d_cost, Nm);
      cudaMalloc((void **)&d_head, Nm);
      cudaMalloc((void **)&d_tail, Nm);

      // Copy once for all
      CHECK(cudaMemcpyAsync(d_cost, G.cost, Nm, cudaMemcpyHostToDevice));
      CHECK(cudaMemcpyAsync(d_head, G.head, Nm, cudaMemcpyHostToDevice));
      CHECK(cudaMemcpyAsync(d_tail, G.tail, Nm, cudaMemcpyHostToDevice));

      // Size for dual variables
      size_t Npi = N * sizeof(int);
      int *h_PI;
      int *d_PI;
      CHECK(cudaMallocHost((void **)&h_PI, Npi));
      CHECK(cudaMalloc((void **)&d_PI, Npi));

      // Negative reduced cost variables
      int threads_per_block = 512;
      int number_of_blocks =
          ((G.m + threads_per_block - 1) / threads_per_block) / 8;

      const int BB = 16;

      int subblock = G.m / BB;
      int B0 = 0, B1 = 0, B2 = 0;
      B2 = std::min(B1 + subblock, G.m);
      size_t Bvar = number_of_blocks / BB * sizeof(int);

      // fprintf(stdout, "blocks=%d, threadsXblock=%d, subblocks=%d ",
      //        number_of_blocks, threads_per_block, number_of_blocks / BB);
      // fflush(stdout);

      size_t Nvar = number_of_blocks / BB * sizeof(int);
      int *h_Var;
      int *d_Var;
      CHECK(cudaMallocHost((void **)&h_Var, Nvar));
      CHECK(cudaMalloc((void **)&d_Var, Nvar));

      Vars vnew;
      vnew.reserve(number_of_blocks);

#// Build the graph for min cost flow
      NetSimplex simplex('E', G.n, 0);

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      int viol = 0;
      simplex.setOptTolerance(-viol);

      // add first d source nodes
      for (int i = 0; i < G.n; ++i) simplex.addNode(i, G.getSupply(i));

      // Init the simplex
      _status = simplex.run();

      while (_status != ProblemType::TIMELIMIT) {
        // Take the dual values
        for (int j = 0; j < N; ++j) h_PI[j] = -simplex.potential(j);

        CHECK(cudaMemcpy(d_PI, h_PI, Npi, cudaMemcpyHostToDevice));

        B0 = B1;
        int itit = 0;
        while (true) {
          deviceReduceKernel<<<number_of_blocks / BB, threads_per_block>>>(
              B1, B2, d_cost, d_head, d_tail, d_PI, d_Var, -viol);

          CHECK(cudaMemcpy(h_Var, d_Var, Bvar, cudaMemcpyDeviceToHost));

          // Take all negative reduced cost variables
          vnew.clear();
          for (int h = 0, h_max = number_of_blocks / BB; h < h_max; ++h)
            if (h_Var[h] > -1)
              vnew.emplace_back(G.head[h_Var[h]], G.tail[h_Var[h]],
                                G.cost[h_Var[h]]);

          if (vnew.empty()) {
            B1 = (B1 + subblock);
            if (B1 == G.m) B1 = 0;
            if (B1 == B0) break;
            B2 = std::min(B1 + subblock, G.m);
          } else {
            B1 = (B1 + subblock);
            if (B1 == G.m) B1 = 0;
            B2 = std::min(B1 + subblock, G.m);
            break;
          }
        }
        if (vnew.empty()) break;

        int new_arcs = simplex.updateArcs(vnew);
        _status = simplex.reRun();

        // if (vnew.empty()) {
        //  if (viol == 0) {

        //  }
        //  if (viol > 1) {
        //    viol = viol >> 1;
        //    // double tt =
        //    //
        //    double(std::chrono::duration_cast<std::chrono::milliseconds>(
        //    //               std::chrono::steady_clock::now() - start_t)
        //    //               .count()) /
        //    //    1000;
        //    // fprintf(stdout, "min viol: %d, cost: %f, time: %f\n", viol,
        //    //        (double)simplex.totalCost() / double(G.n), tt);
        //  } else
        //    viol = 0;
        //  simplex.setOptTolerance(-viol);
        //}

        // std::sort(vnew.begin(), vnew.end(),
        //          [](const Var &v, const Var &w) { return v.c > w.c; });

        // Replace old constraints with new ones
        // int new_arcs = simplex.updateArcs(vnew);

        //_status = simplex.reRun();
      }

      double fobj = (double)simplex.totalCost() / double(G.n);

      _iterations = (int)simplex.iterations();
      _runtime = simplex.runtime();
      _num_arcs = simplex.num_arcs();

      auto end_t = std::chrono::steady_clock::now();
      auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end_t - start_t)
                             .count()) /
                  1000;

      PRINT(
          "CUDALGEN %s it %d UB %.6f runtime %.4f simplex %.4f "
          "num_arcs %d blocks %d\n",
          msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs,
          number_of_blocks);

      cudaFree(d_cost);
      cudaFree(d_head);
      cudaFree(d_tail);
      cudaFree(h_PI);
      cudaFree(d_PI);
      cudaFree(h_Var);
      cudaFree(d_Var);

      return _all;
#else
      fprintf(stdout, "Kernel CUDA not implemented.\n");

      return -1;
#endif
    }
    //#endif

    if (filename == "colgen") {
      auto start_t = std::chrono::steady_clock::now();

      int N = G.n;

      // Dual multipliers
      vector<int> pi;
      pi.resize(N, 0);

      // Auxiliary variables
      int NUM_BLOCKS = 2 * N;
      Vars vars(NUM_BLOCKS);  // one for block o threads?

      // TODO: pensare a come gestire questo caso!
      // for (int i = 0; i < G.n; ++i)
      //  vars[i].a = i;

      Vars vnew;
      vnew.reserve(G.n);

      // Build the graph for min cost flow
      NetSimplex simplex('E', G.n, 0);

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      simplex.setOptTolerance(0);

      // add first d source nodes
      for (int i = 0; i < G.n; ++i) simplex.addNode(i, G.getSupply(i));

      // Init the simplex
      _status = simplex.run();

      std::vector<int> indexes;
      indexes.resize(G.m);
      std::iota(indexes.begin(), indexes.end(), 0);
      std::random_device rd;
      std::mt19937 g(rd());

      // Solve separation problem:
      int block = G.m / NUM_BLOCKS;

      int cnt = 0;
      while (_status != ProblemType::TIMELIMIT) {
        // Take the dual values
        for (int j = 0; j < N; ++j) pi[j] = -simplex.potential(j);

        cnt++;

        // TODO: I can try a thread pool instead OMP
#pragma omp parallel for
        for (int b = 0; b < NUM_BLOCKS; ++b) {
          int best_e = 0;
          int best_v = 0;
          int best_c = std::numeric_limits<int>::max();

          for (int e = b * block, e_max = std::min<int>(G.m, (b + 1) * block);
               e < e_max; ++e) {
            int violation = G.cost[e] - pi[G.head[e]] + pi[G.tail[e]];
            if (violation < best_v ||
                (violation == best_v && G.cost[e] < best_c)) {
              best_e = e;
              best_v = violation;
              best_c = G.cost[e];
              if (best_v < 0) break;
            }
          }

          // Store most violated cuts for element i
          if (best_v < 0) {
            vars[b].a = G.head[best_e];
            vars[b].b = G.tail[best_e];
            vars[b].c = G.cost[best_e];
          } else
            vars[b].c = -1;
        }

        // Take all negative reduced cost variables
        vnew.clear();
        for (int i = 0; i < NUM_BLOCKS; i++) {
          if (vars[i].c > -1) vnew.push_back(vars[i]);
          vars[i].c = -1;
        }

        if (vnew.empty()) break;

        // std::sort(vnew.begin(), vnew.end(),
        //          [](const Var &v, const Var &w) { return v.c < w.c; });

        // Replace old constraints with new ones
        int new_arcs = simplex.updateArcs(vnew);

        _status = simplex.reRun();
      }

      double fobj = (double)simplex.totalCost() / double(G.n);

      _iterations = (int)simplex.iterations();
      _runtime = simplex.runtime();
      _num_arcs = simplex.num_arcs();

      auto end_t = std::chrono::steady_clock::now();
      auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end_t - start_t)
                             .count()) /
                  1000;

      PRINT(
          "MYCOLGEN %s it %d UB %.6f runtime %.4f simplex %.4f "
          "num_arcs %d\n",
          msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

      return _all;
    }

    if (filename == "bipartiteEati") {
      auto start_t = std::chrono::steady_clock::now();

      int N = G.n;

      // Build the graph for min cost flow
      NetSimplex simplex('F', G.n, G.m);

      // Set the parameters
      simplex.setTimelimit(timelimit);
      simplex.setVerbosity(verbosity);
      simplex.setOptTolerance(0);

      // add first d source nodes
      for (int i = 0; i < G.n; ++i) simplex.addNode(i, G.getSupply(i));

      for (int e = 0; e < G.m; ++e)
        simplex.addArc(G.head[e], G.tail[e], G.cost[e]);

      // Init the simplex
      simplex.run();

      double fobj = (double)simplex.totalCost() / double(G.n);

      _iterations = (int)simplex.iterations();
      _runtime = simplex.runtime();
      _num_arcs = simplex.num_arcs();

      auto end_t = std::chrono::steady_clock::now();
      auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end_t - start_t)
                             .count()) /
                  1000;

      PRINT(
          "BIP-EATI %s it %d UB %.6f runtime %.4f simplex %.4f "
          "num_arcs %d\n",
          msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

      return _all;
    }

    if (filename == "bipartiteEapi") {
      auto start_t = std::chrono::steady_clock::now();

      int N = G.n;

      // Build the graph for min cost flow
      EapiSimplex simplex(G.n, G.m);

      // add first d source nodes
      for (int i = 0; i < G.n; ++i) simplex.addNode(i, G.getSupply(i));

      for (int e = 0; e < G.m; ++e)
        simplex.addArc(G.head[e], G.tail[e], G.cost[e]);

      // Init the simplex
      double fobj = simplex.solve();

      fobj = fobj / double(G.n);

      _iterations = (int)simplex.iterations();
      _runtime = simplex.runtime();
      _num_arcs = simplex.M - simplex.M0;

      auto end_t = std::chrono::steady_clock::now();
      auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                             end_t - start_t)
                             .count()) /
                  1000;

      PRINT(
          "BIP-EAPI %s it %d UB %.6f runtime %.4f simplex %.4f "
          "num_arcs %d\n",
          msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

      return _all;
    }

    //--------------------------------------------
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
      CPXsetdblparam(env, CPXPARAM_Network_Tolerances_Feasibility, 1e-09);
      CPXsetdblparam(env, CPXPARAM_Network_Tolerances_Optimality, 1e-09);

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
      for (int i = 0; i < G.n; ++i) ss[i] = -G.getSupply(i);

      CPXNETaddnodes(env, net, (int)G.n, &ss[0], NULL);

      // Add all arcs
      vector<int> tail;
      vector<int> head;
      vector<double> obj(G.m, 0);
      for (int e = 0; e < G.m; ++e) obj[e] = G.cost[e];

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

      if (solstat != CPX_STAT_OPTIMAL)
        fprintf(stdout, "ERROR: Cplex status: %d\n", solstat);

      double objval = 0.0;
      CPXNETgetobjval(env, net, &objval);

      objval = objval / double(G.n);

      int _iterations = CPXNETgetitcnt(env, net);

      /* Free up the problem as allocated by CPXNETcreateprob, if necessary
       */
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

      PRINT(
          "BIP-PLEX %s it %d UB %.6f runtime %.4f simplex %.4f "
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
};

}  // namespace DOT
