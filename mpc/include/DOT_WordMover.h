/**
 * @fileoverview Copyright (c) 2019-2023, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include "DOT_Commons.h"
#include "DOT_Cplex.h"
#include "DOT_NetSimplexDouble.h"

// NAME SPACE
namespace DOT {

class WordMover {
 public:
  // Standard c'tor
  WordMover() : n(0), m(0), q(0) {}

  // copy c'tor
  WordMover(const WordMover &o) : n(o.n), m(o.m), q(o.q) {
    X = (double *)malloc(n * q * sizeof(double));
    H = (double *)malloc(m * n * sizeof(double));
  }

  // desctuctor
  ~WordMover() {
    if (X != nullptr) free(X);
    if (H != nullptr) free(H);
    if (C != nullptr) free(C);
  }

  // parse embedding vector from data from file
  void parseX(const std::string &filename) {
    std::ifstream infile(filename.c_str(),
                         std::ios::binary | std::ios::in | std::ios::ate);
    if (!infile.is_open()) {
      fprintf(stdout, "Error opening file!\n");
      exit(0);
    }

    size_t size = infile.tellg();

    char *memblock = (char *)malloc(size * sizeof(char));
    infile.seekg(0, std::ios::beg);
    infile.read(memblock, size);
    infile.close();

    double *double_values = (double *)memblock;  // reinterpret as doubles
    n = (size_t)double_values[0];                // Number of vectors
    q = (size_t)double_values[1];                // Size of vectors
    size_t Xsize = n * q;

    X = (double *)malloc(Xsize * sizeof(double));
    for (size_t i = 0; i < Xsize; ++i) X[i] = (double)double_values[2 + i];

    free(memblock);
  }

  // parse embedding vector from data from file
  void parseH(const std::string &filename) {
    std::ifstream infile(filename.c_str(),
                         std::ios::binary | std::ios::in | std::ios::ate);
    if (!infile.is_open()) {
      fprintf(stdout, "Error opening file!\n");
      exit(0);
    }

    size_t size = infile.tellg();

    char *memblock = (char *)malloc(size * sizeof(char));
    infile.seekg(0, std::ios::beg);
    infile.read(memblock, size);
    infile.close();

    double *double_values = (double *)memblock;  // reinterpret as doubles
    m = (size_t)double_values[0];                // number of text documents
    n = (size_t)double_values[1];  // size of vectors support the documents
    size_t Hsize = m * n;

    H = (double *)malloc(Hsize * sizeof(double));
    for (size_t i = 0; i < Hsize; ++i) H[i] = (double)double_values[2 + i];

    free(memblock);
  }

  // Compute distances between matrices
  void computeCostMatrix() {
    C = (double *)malloc(n * n * sizeof(double));
    for (size_t i = 0; i < n; ++i) {
      C[i * n + i] = 0.0;
      for (size_t j = i + 1; j < n; ++j) {
        size_t idx = i * n + j;
        C[idx] = 0.0;
        for (size_t h = 0; h < q; ++h)
          C[idx] += pow(X[i * q + h] - X[j * q + h], 2);
        C[idx] = sqrt(C[idx]);
        C[j * n + i] = C[idx];
      }
    }
  }

  // Normalize histograms
  void normalizeHistogramRow(size_t i) {
    double s = 0.0;
    for (size_t h = 0; h < n; h++) s += H[i * n + h];
    for (size_t h = 0; h < n; h++) H[i * n + h] = H[i * n + h] / s;
  }

  void normalizeAllHistograms(void) {
    for (size_t i = 0; i < m; ++i) normalizeHistogramRow(i);
  }

  // Compute distance between a pair of text documents
  double distanceNetsimplex(size_t a, size_t b, const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    vector<int> IDhead;
    vector<int> IDtail;
    IDhead.reserve(n);
    IDtail.reserve(n);

    vector<double> supply;
    supply.reserve(2 * n);
    for (int i = 0; i < n; ++i)
      if (H[a * n + i] > 0.0) {
        IDtail.push_back(i);
        supply.push_back(H[a * n + i]);
      }

    int N = supply.size();
    fprintf(stdout, "H1 size: %d\n", (int)N);

    for (int j = 0; j < n; ++j)
      if (H[b * n + j] > 0.0) {
        IDhead.push_back(j);
        supply.push_back(-H[b * n + j]);  // Destinations are negative values
      }
    fprintf(stdout, "H2 size: %d\n", (int)supply.size() - N);

    int Na = IDtail.size();
    int Nb = IDhead.size();
    int Nab = Na * Nb;

    // Build the graph for min cost flow
    NetSimplexDouble simplex('F', Na + Nb, Na * Nb);

    for (int i = 0, i_max = Na + Nb; i < i_max; ++i) {
      simplex.addNode(i, supply[i]);
      fprintf(stdout, "n %d %f\n", i, supply[i]);
    }

    for (int i = 0; i < Na; ++i)
      for (int j = 0; j < Nb; ++j) {
        simplex.addArc(i, Na + j, C[IDtail[i] * n + IDhead[j]]);
        fprintf(stdout, "%d %d %f\n", i, Na + j, C[IDtail[i] * n + IDhead[j]]);
      }

    // Set the parameters
    // simplex.setTimelimit(timelimit);
    // simplex.setVerbosity(verbosity);
    // simplex.setOptTolerance(opt_tolerance);

    auto _status = simplex.run();

    double _runtime = simplex.runtime();
    uint64_t _iterations = simplex.iterations();
    uint64_t _num_arcs = simplex.num_arcs();
    uint64_t _num_nodes = simplex.num_nodes();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    double fobj = simplex.totalCost();

    // Upper bound on the missed mass
    DOT::logger.info(
        "MYCOLGEN %s it %d LB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return fobj;
  }

  // Compute distance between a pair of text documents
  double distanceColumnGeneration(size_t a, size_t b,
                                  const std::string &msg = "") {
    auto start_t = std::chrono::steady_clock::now();

    vector<int> IDhead;
    vector<int> IDtail;
    IDhead.reserve(n);
    IDtail.reserve(n);
    vector<int> INVhead(n, 0);
    vector<int> INVtail(n, 0);
    vector<double> supply;
    supply.reserve(2 * n);
    for (int i = 0; i < n; ++i)
      if (H[a * n + i] > 0.0) {
        IDhead.push_back(supply.size());
        INVhead[supply.size()] = i;
        supply.push_back(H[a * n + i]);
      }

    size_t N = supply.size();
    fprintf(stdout, "H1 size: %d\n", (int)N);

    for (int j = 0; j < n; ++j)
      if (H[b * n + j] > 0.0) {
        IDtail.push_back(supply.size());
        INVtail[supply.size()] = j;
        supply.push_back(-H[b * n + j]);  // Destinations are negative values
      }
    fprintf(stdout, "H2 size: %d\n", (int)supply.size() - N);

    size_t Na = IDhead.size();
    size_t Nb = IDtail.size();
    size_t Nab = Na * Nb;

    // STEGUA: E' sbagliato Na!!!
    vector<double> pi(Na + Nb, 0.0);
    pi.shrink_to_fit();

    DoubleVars vars(Na);
    for (size_t i = 0; i < Na; ++i) vars[i].a = i;

    DoubleVars vnew;
    vnew.reserve(Na);

    // Build the graph for min cost flow
    NetSimplexDouble simplex('E', Na + Nb, 0);

    for (size_t i = 0, i_max = Na + Nb; i < i_max; ++i)
      simplex.addNode(i, supply[i]);

    // Set the parameters
    // simplex.setTimelimit(timelimit);
    // simplex.setVerbosity(verbosity);
    // simplex.setOptTolerance(opt_tolerance);

    auto _status = simplex.run();

    while (_status != ProblemType::TIMELIMIT) {
      // Take the dual values
      for (int j = 0; j < Na + Nb; ++j) pi[j] = -simplex.potential(j);

      // Solve separation problem
      //#pragma omp parallel
      for (size_t i = Na; i < Na; ++i) {
        double best_v = -1e-09;
        double best_c = -10.0;
        int best_j = 0;  // best second node

        for (size_t j = 0; j < Nb; ++j) {
          double violation = C[INVhead[i] * n + INVtail[j]] - pi[i] + pi[N + j];
          if (violation < best_v) {
            best_v = violation;
            best_j = N + j;
            best_c = C[INVhead[i] * n + INVtail[j]];
          }
        }

        // Store most violated cuts for element i
        vars[i].b = best_j;
        vars[i].c = best_c;
      }

      // Take all negative reduced cost variables
      vnew.clear();
      for (auto &v : vars) {
        if (v.c > -1.0) vnew.push_back(v);
        v.c = -10.0;
      }

      if (vnew.empty()) {
        fprintf(stdout, "no new arcs\n");
        break;
      }

      std::sort(
          vnew.begin(), vnew.end(),
          [](const DoubleVar &v, const DoubleVar &w) { return v.c > w.c; });

      // Replace old constraints with new ones
      int new_arcs = simplex.updateArcs(vnew);

      _status = simplex.reRun();
    }

    double _runtime = simplex.runtime();
    uint64_t _iterations = simplex.iterations();
    uint64_t _num_arcs = simplex.num_arcs();
    uint64_t _num_nodes = simplex.num_nodes();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    double fobj = simplex.totalCost();

    // Upper bound on the missed mass
    DOT::logger.info(
        "MYCOLGEN %s it %d LB %.6f runtime %.4f simplex %.4f "
        "num_arcs %d",
        msg.c_str(), (int)_iterations, fobj, _all, _runtime, (int)_num_arcs);

    return fobj;
  }

  // Compute distance between a pair of text documents
  double distanceCplex(size_t a, size_t b, const std::string &msg = "") {
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

    CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
    CPXsetdblparam(env, CPXPARAM_Network_Tolerances_Feasibility, 1e-09);
    CPXsetdblparam(env, CPXPARAM_Network_Tolerances_Optimality, 1e-09);

    net = CPXNETcreateprob(env, &status, "netCplex");

    if (net == NULL) {
      fprintf(stderr, "Failed to create network object.\n");
      fflush(stderr);
      exit(-1);
    }

    // Add all arcs
    vector<int> IDhead(n, 0);
    vector<int> IDtail(n, 0);
    vector<double> supply;
    supply.reserve(2 * n);
    for (int i = 0; i < n; ++i)
      if (H[a * n + i] > 0.0) {
        IDhead[i] = supply.size();
        supply.push_back(H[a * n + i]);
      }

    size_t N = supply.size();
    fprintf(stdout, "H1 size: %d\n", (int)N);

    for (int j = 0; j < n; ++j)
      if (H[b * n + j] > 0.0) {
        IDtail[j] = supply.size();
        supply.push_back(-H[b * n + j]);  // Destinations are negative values
      }
    fprintf(stdout, "H2 size: %d\n", (int)supply.size() - N);

    vector<int> tail;
    vector<int> head;
    vector<double> obj;

    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        if (H[a * n + i] > 0.0 && H[b * n + j] > 0.0) {
          tail.push_back(IDhead[i]);
          head.push_back(IDtail[j]);
          obj.push_back(C[i * n + j]);
          fprintf(stdout, "%d %d %f\n", IDhead[i], IDtail[j], C[i * n + j]);
        }

    // Build Network instance
    CPXNETchgobjsen(env, net, CPX_MIN);

    CPXNETaddnodes(env, net, (int)supply.size(), &supply[0], NULL);

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

    fprintf(stdout, "status: %d, solstat: %d\n", status, solstat);

    double objval = 0.0;
    CPXNETgetobjval(env, net, &objval);

    int _iterations = CPXNETgetitcnt(env, net);

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

    fprintf(stdout, "BIP-PLEX %s it %d UB %.6f runtime %.4f simplex %.4f ",
            msg.c_str(), (int)_iterations, objval, _all, _runtime);
  }

  // Dump the dimensions
  void dump() {
    fprintf(stdout, "Dimensions -> m: %d, n: %d, q: %d\n", (int)m, (int)n,
            (int)q);
  }
  void dumpRowX(size_t i) {
    fprintf(stdout, "X[%d] = ", i);
    for (size_t h = 0; h < q; ++h) fprintf(stdout, "%f ", X[i * q + h]);
    fprintf(stdout, "\n");
  }

  void dumpRowH(size_t i) {
    fprintf(stdout, "H[%d] = ", i);
    for (size_t h = 0; h < n; ++h) fprintf(stdout, "%f ", H[i * n + h]);
    fprintf(stdout, "\n");
  }

 private:
  size_t m;  // number of text document
  size_t n;  // number of vector x_i in the embedding X
  size_t q;  // Size of the embedding space x_i \in R^q

  double *X;  // Embedding vectors
  double *H;  // Histograms of text document
  double *C;  // Cost matrix of dimension n*n
};

}  // namespace DOT