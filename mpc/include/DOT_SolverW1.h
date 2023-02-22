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

//#include "yocta_logger.hh"

namespace std {
typedef tuple<int, int, int, int> tuple4int;

typedef unordered_map<tuple4int, int, hash_tuple4int> var2idx;
}  // namespace std

typedef std::pair<int, int> int_pair;

// Zeta block for coordinates vector
#define BLOCKSIZE 1024
#define BLOCKNUM 8

namespace DOT {

// Solver class, which wrapper the Network Simplex algorithm
class SolverW1 {
 public:
  // Standard c'tor
  SolverW1()
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
  double cplex(const Histogram2D &A, const Histogram2D &B,
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

    vector<double> supply(n * n, 0);

    auto ID = [&n](int x, int y) { return x * n + y; };

    // Add all nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) supply[ID(i, j)] = A.get(i, j) - B.get(i, j);

    CPXNETaddnodes(env, net, n * n, &supply[0], NULL);

    // Add all arcs
    vector<int> tail;
    tail.reserve(20 * n * n);
    vector<int> head;
    head.reserve(20 * n * n);
    vector<double> obj;
    obj.reserve(20 * n * n);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
            tail.push_back(ID(i, j));
            head.push_back(ID(i + v, j + w));
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
        "FLOCPLEX %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, objval, _all, _runtime, (int)_num_arcs);

    return _all;
  }

  //--------------------------------------------------------------------------
  std::array<double, 4> eapi(const Histogram2D &A, const Histogram2D &B,
                             const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    int mm = 0;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) mm++;
        }
      }

    // Build the graph for min cost flow
    EapiSimplex simplex(n * n, n * n + mm);

    auto ID = [&n](int x, int y) { return x * n + y; };

    // Add all nodes
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

    DOT::logger.info(
        "FLO-EAPI %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)mm);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
  }

  //--------------------------------------------------------------------------
  std::array<double, 4> eati(const Histogram2D &A, const Histogram2D &B,
                             const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    // Build the graph for min cost flow
    NetSimplex simplex('F', static_cast<int>(n * n),
                       static_cast<int>(n * n) * static_cast<int>(n * n));

    // Set the parameters
    simplex.setTimelimit(timelimit);
    simplex.setVerbosity(verbosity);
    simplex.setOptTolerance(0);

    auto ID = [&n](int x, int y) { return x * n + y; };

    // Add all nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(ID(i, j), A.get(i, j) - B.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n)
            simplex.addArc(ID(i, j), ID(i + v, j + w), p.c_vw);
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

    DOT::logger.info(
        "FLO-EATI %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
    ;
  }

  double gurobi(const Histogram2D &A, const Histogram2D &B,
                const std::string &msg = "", int method = 1,
                bool dump = false) {
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

    GRBsetdblparam(env, GRB_DBL_PAR_FEASIBILITYTOL, 1e-09);
    GRBsetdblparam(env, GRB_DBL_PAR_OPTIMALITYTOL, 1e-09);

    auto ID = [&n](int x, int y) { return x * n + y; };

    vector<double> obj;

    // First layer move horizontal
    std::var2idx X;
    int idx = 0;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
            obj.push_back(p.c_vw);
            std::tuple4int p = std::tuple4int(i, j, i + v, j + w);
            X[p] = idx;
            idx++;
          }
        }
      }

    int m = (int)obj.size();

    // Variables and objective function
    GRBnewmodel(env, &model, "ot", m, &obj[0], NULL, NULL, NULL, NULL);

    // Add constraints
    vector<double> val(n, 1.0);
    vector<int> ind;
    ind.reserve(n);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        ind.clear();
        val.clear();
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
            ind.push_back(X.at(std::tuple(i, j, i + v, j + w)));
            val.push_back(1.0);
            ind.push_back(X.at(std::tuple(i + v, j + w, i, j)));
            val.push_back(-1.0);
          }
        }

        GRBaddconstr(model, (int)ind.size(), &ind[0], &val[0], GRB_EQUAL,
                     A.get(i, j) - B.get(i, j), NULL);
      }

    if (dump) GRBwrite(model, (msg + "rlp.gz").c_str());

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
        "FLO-GUR%d %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        method, msg.c_str(), (int)_iterations, objval, _all, _runtime,
        (int)_num_arcs);

    return _runtime;
  }

  double lemon(const Histogram2D &A, const Histogram2D &B,
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
    // Add all nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) nodes_A.emplace_back(g.addNode());

    // Add arcs for complete bipartite graph
    std::vector<Graph::Arc> arcs;
    arcs.reserve(20 * d);
    std::vector<int64_t> arcs_costs;
    arcs_costs.reserve(20 * d);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
            arcs.emplace_back(
                g.addArc(nodes_A[ID(i, j)], nodes_A[ID(i + v, j + w)]));
            arcs_costs.emplace_back(p.c_vw);
          }
        }
      }

    // lower and upper bounds, cost
    ListDigraph::ArcMap<int64_t> l_i(g), u_i(g);
    ListDigraph::ArcMap<int64_t> c_i(g);

    // FLow balance
    ListDigraph::NodeMap<int64_t> b_i(g);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        b_i[nodes_A[ID(i, j)]] = A.get(i, j) - B.get(i, j);
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

    auto ss_t = std::chrono::steady_clock::now();
    NetworkSimplex<Graph, int64_t, int64_t>::ProblemType ret = simplex.run();
    auto ee_t = std::chrono::steady_clock::now();
    _runtime = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                          ee_t - ss_t)
                          .count()) /
               1000;

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
        "FLO-LEMO %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return _runtime;
  }

  //--------------------------------------------------------------------------
  void init_coprimes(int L, int N) {
    // La distanza la normalizzo in modo che la distanza massima sia
    // paria  1.0 e poi la riscalo all'intero più vicino dopo aver
    // moltiplicato per un fattore di scala

    double Dmax = 2 * L * L;

    coprimes.clear();
    for (int v = -L; v <= L; ++v)
      for (int w = -L; w <= L; ++w)
        if (!(v == 0 && w == 0) && std::gcd(v, w) == 1)
          coprimes.emplace_back(
              v, w, int(round(SCALE * sqrt(double(v * v + w * w) / Dmax))));
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
