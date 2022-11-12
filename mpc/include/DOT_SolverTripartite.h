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

namespace std {
typedef tuple<int, int, int, int> tuple4int;

struct hash_tuple4int {
  inline size_t operator()(const tuple4int &x) const {
    std::hash<int> int_hasher;
    return int_hasher(std::get<0>(x)) ^ int_hasher(std::get<1>(x)) ^
           int_hasher(std::get<2>(x)) ^ int_hasher(std::get<3>(x));
  }
};

typedef unordered_map<tuple4int, int, hash_tuple4int> var2idx;
}  // namespace std

namespace DOT {

// Solver class, which wrapper the Network Simplex algorithm
class SolverTripartite {
 public:
  // Standard c'tor
  SolverTripartite()
      : _runtime(0.0),
        _n_log(0),
        L(-1),
        verbosity(DOT_VAL_INFO),
        recode(""),
        opt_tolerance(0),
        timelimit(std::numeric_limits<double>::max()) {}

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

  //--------------------------------------------------------------------------
  std::array<double, 4> eapi(const Histogram2D &A, const Histogram2D &B,
                             const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    // Build the graph for min cost flow
    EapiSimplex simplex(static_cast<int>(3 * n * n),
                        static_cast<int>(2 * n * n * n));

    auto ID = [&n](int x, int y) { return x * n + y; };

    // add nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) simplex.addNode(ID(i, j), A.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) simplex.addNode(n * n + ID(i, j), 0.0);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(2 * n * n + ID(i, j), -B.get(i, j));

    // add arcs
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        for (int v = 0; v < n; ++v)
          simplex.addArc(ID(i, j), n * n + ID(i, v), pow(j - v, 2));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        for (int v = 0; v < n; ++v)
          simplex.addArc(n * n + ID(i, j), 2 * n * n + ID(v, j), pow(i - v, 2));

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
        "TRI-EAPI %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
  }

  //--------------------------------------------------------------------------
  std::array<double, 4> eati(const Histogram2D &A, const Histogram2D &B,
                             const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    // Build the graph for min cost flow
    NetSimplex simplex('F', static_cast<int>(3 * n * n),
                       static_cast<int>(2 * n * n * n));

    // Set the parameters
    simplex.setTimelimit(timelimit);
    simplex.setVerbosity(verbosity);
    simplex.setOptTolerance(0);

    auto ID = [&n](int x, int y) { return x * n + y; };

    // add nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) simplex.addNode(ID(i, j), A.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) simplex.addNode(n * n + ID(i, j), 0.0);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(2 * n * n + ID(i, j), -B.get(i, j));

    // add arcs
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        for (int v = 0; v < n; ++v)
          simplex.addArc(ID(i, j), n * n + ID(i, v), pow(j - v, 2));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        for (int v = 0; v < n; ++v)
          simplex.addArc(n * n + ID(i, j), 2 * n * n + ID(v, j), pow(i - v, 2));

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
        "TRI-EATI %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return {_runtime, simplex._time_pricing, simplex._time_update_basis,
            simplex._time_update_duals};
    ;
  }

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

    vector<double> supply(3 * n * n, 0);

    auto ID = [&n](int x, int y) { return x * n + y; };

    // Add all nodes
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) supply[ID(i, j)] = A.get(i, j);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) supply[n * n + ID(i, j)] = 0.0;

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) supply[2 * n * n + ID(i, j)] = -B.get(i, j);

    CPXNETaddnodes(env, net, 3 * n * n, &supply[0], NULL);

    // Add all arcs
    vector<int> tail;
    tail.reserve(2 * n * n * n);
    vector<int> head;
    head.reserve(2 * n * n * n);
    vector<double> obj;
    obj.reserve(2 * n * n * n);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        for (int v = 0; v < n; ++v) {
          tail.push_back(ID(i, j));
          head.push_back(n * n + ID(i, v));
          obj.push_back(pow(j - v, 2));
        }
      }

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        for (int v = 0; v < n; ++v) {
          tail.push_back(n * n + ID(i, j));
          head.push_back(2 * n * n + ID(v, j));
          obj.push_back(pow(i - v, 2));
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
        "TRI-PLEX %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, objval, _all, _runtime, (int)_num_arcs);

    return _all;
  }

  double gurobi(const Histogram2D &A, const Histogram2D &B,
                const std::string &msg = "", int method = 3) {
    // Init coprimes data structure
    int n = A.getN();

    auto start_t = std::chrono::steady_clock::now();

    GRBenv *env = NULL;
    GRBmodel *model = NULL;
    int status = 0;

    GRBloadenv(&env, NULL);
    GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);
    // GRBsetintparam(env, GRB _INT_PAR_THREADS, 1);
    if (method == 1) GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
    if (method == 2) {
      GRBsetintparam(env, GRB_INT_PAR_CROSSOVER, 0);
      GRBsetintparam(env, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER);
    }
    if (method == 3) GRBsetintparam(env, GRB_INT_PAR_NETWORKALG, 1);

    GRBsetdblparam(env, GRB_DBL_PAR_FEASIBILITYTOL, 1e-09);
    GRBsetdblparam(env, GRB_DBL_PAR_OPTIMALITYTOL, 1e-09);

    vector<double> obj;
    vector<double> lb;

    // First layer move horizontal
    std::var2idx X;
    std::var2idx Y;
    int idx = 0;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        // From (i,j) to N+(v, j)
        for (int v = 0; v < n; ++v) {
          obj.push_back(pow(j - v, 2));
          std::tuple4int p = std::tuple4int(i, j, i, v);
          X[p] = idx;
          idx++;
        }
      }

    // Second layer move vertically
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        // From N+(i,j) to 2*N+(i, w)
        for (int v = 0; v < n; ++v) {
          obj.push_back(pow(i - v, 2));
          std::tuple4int p = std::tuple4int(i, j, v, j);
          Y[p] = idx;
          idx++;
        }
      }

    obj.shrink_to_fit();
    int m = (int)obj.size();

    // Variables and objective function
    GRBnewmodel(env, &model, "ot", m, &obj[0], NULL, NULL, NULL, NULL);

    // Add constraints
    vector<double> val(2 * n, 1.0);
    vector<int> ind;
    ind.reserve(2 * n);

    // print('2. Add initial constraint sets')
    // for i,j in P:
    //    m.addConstr(quicksum(X[i,j,i,v] for v in range(n)) <= h1[i,j])
    // for i,j in P:
    //    m.addConstr(quicksum(Y[v,j,i,j] for v in range(n)) >= h2[i,j])
    // for i,j in P:
    //    m.addConstr(quicksum(X[i,v,i,j] for v in range(n))-quicksum(Y[i,j,v,j]
    //    for v in range(n)) == 0)

    // First Layer
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        ind.clear();
        for (int v = 0; v < n; ++v)
          ind.push_back(X.at(std::tuple4int(i, j, i, v)));

        GRBaddconstr(model, n, &ind[0], &val[0], GRB_EQUAL, A.get(i, j), NULL);
      }

    // Middle layer
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        ind.clear();
        for (int v = 0; v < n; ++v) {
          ind.push_back(X.at(std::tuple4int(i, v, i, j)));
          val[v] = 1.0;
        }
        for (int v = 0; v < n; ++v) {
          ind.push_back(Y.at(std::tuple4int(i, j, v, j)));
          val[n + v] = -1.0;
        }

        GRBaddconstr(model, 2 * n, &ind[0], &val[0], GRB_EQUAL, 0, NULL);
      }

    // Last layer
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        ind.clear();
        for (int v = 0; v < n; ++v)
          ind.push_back(Y.at(std::tuple4int(v, j, i, j)));

        GRBaddconstr(model, n, &ind[0], &val[0], GRB_EQUAL, B.get(i, j), NULL);
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
        "TRI-GUR%d %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        method, msg.c_str(), (int)_iterations, objval, _all, _runtime,
        (int)_num_arcs);

    return _all;
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
    // add first d source nodes for first partition
    for (size_t i = 0; i < d; ++i) nodes_A.emplace_back(g.addNode());

    std::vector<Graph::Node> nodes_B;
    nodes_B.reserve(d);
    // add first d source nodes for second partition
    for (size_t i = 0; i < d; ++i) nodes_B.emplace_back(g.addNode());

    std::vector<Graph::Node> nodes_C;
    nodes_C.reserve(d);
    // add first d source nodes for second partition
    for (size_t i = 0; i < d; ++i) nodes_C.emplace_back(g.addNode());

    // Add arcs for complete bipartite graph
    std::vector<Graph::Arc> arcs;
    arcs.reserve(2 * d * n);
    std::vector<int64_t> arcs_costs;
    arcs_costs.reserve(2 * d * n);

    // add arcs
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        for (int v = 0; v < n; ++v) {
          arcs.emplace_back(g.addArc(nodes_A[ID(i, j)], nodes_B[ID(i, v)]));
          arcs_costs.emplace_back(pow(j - v, 2));
        }

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        for (int v = 0; v < n; ++v) {
          arcs.emplace_back(g.addArc(nodes_B[ID(i, j)], nodes_C[ID(v, j)]));
          arcs_costs.emplace_back(pow(i - v, 2));
        }

    // lower and upper bounds, cost
    ListDigraph::ArcMap<int64_t> l_i(g), u_i(g);
    ListDigraph::ArcMap<int64_t> c_i(g);

    // FLow balance
    ListDigraph::NodeMap<int64_t> b_i(g);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) {
        b_i[nodes_A[ID(i, j)]] = +A.get(i, j);
        b_i[nodes_B[ID(i, j)]] = 0.0;
        b_i[nodes_C[ID(i, j)]] = -B.get(i, j);
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
        "TRI-LEMO %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return _all;
  }

  double colgen(const Histogram2D &A, const Histogram2D &B,
                const std::string &msg) {
    auto start_t = std::chrono::steady_clock::now();

    int n = A.getN();

    auto ID = [&n](int x, int y) { return x * n + y; };

    int N = 3 * n * n;
    vector<int> pi(N, 0);
    pi.shrink_to_fit();

    Vars vars(n * n);

    Vars vnew;
    vnew.reserve(n * n);

    // Build the graph for min cost flow
    NetSimplex simplex('E', N, 0);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) simplex.addNode(ID(i, j), A.get(i, j));

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) simplex.addNode(n * n + ID(i, j), 0.0);

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        simplex.addNode(2 * n * n + ID(i, j), -B.get(i, j));

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
          int best_a = -1;
          int best_b = 0;  // best second node
          int h = ID(i, j);

          for (int v = 0; v < n; ++v) {
            int c = pow(j - v, 2);
            int violation = c - pi[h] + pi[n * n + ID(i, v)];
            if (violation < best_v) {
              best_v = violation;
              best_c = c;
              best_b = n * n + ID(i, v);
              best_a = h;
            }
          }

          int hh = n * n + ID(i, j);

          for (int v = 0; v < n; ++v) {
            int c = pow(i - v, 2);
            int violation = c - pi[hh] + pi[2 * n * n + ID(v, j)];
            if (violation < best_v) {
              best_v = violation;
              best_c = c;
              best_b = 2 * n * n + ID(v, j);
              best_a = hh;
            }
          }

          // Store most violated cuts for element i
          vars[h].a = best_a;
          vars[h].b = best_b;
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
        "TRCOLGEN %s it %d UB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return _all;
  }

 private:
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
