/**
 * @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <omp.h>

#include <vector>
using std::vector;

#include <set>
using std::set;

#include <array>
using std::array;

#include <unordered_map>
using std::unordered_map;

#include <unordered_set>
using std::unordered_set;

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifndef __GCD_
#define __GCD_
int GCD(int _a, int _b) {
  int a = (_a >= 0 ? _a : -_a);
  int b = (_b >= 0 ? _b : -_b);
  while (b != 0) {
    int t = b;
    b = a % b;
    a = t;
  }
  return a;
}
#endif

void tolower(std::string &data) {
  std::transform(data.begin(), data.end(), data.begin(),
                 [](unsigned char c) { return std::tolower(c); });
}

// List of solver parameters:
//------------------------------

// (exact, approx)
constexpr auto DOT_PAR_METHOD = "Method";
constexpr auto DOT_VAL_EXACT = "exact";
constexpr auto DOT_VAL_APPROX = "approx";

// (bipartite, mincostflow)
constexpr auto DOT_PAR_MODEL = "Model";
constexpr auto DOT_VAL_BIPARTITE = "bipartite";
constexpr auto DOT_VAL_MINCOSTFLOW = "mincostflow";

// (fullmodel, colgen)
constexpr auto DOT_PAR_ALGORITHM = "Algorithm";
constexpr auto DOT_VAL_FULLMODEL = "fullmodel";
constexpr auto DOT_VAL_COLGEN = "colgen";

// ('SILENT', 'INFO', 'DEBUG')
constexpr auto DOT_PAR_VERBOSITY = "Verbosity";
constexpr auto DOT_VAL_SILENT = "silent";
constexpr auto DOT_VAL_INFO = "info";
constexpr auto DOT_VAL_DEBUG = "debug";

constexpr auto DOT_PAR_TIMELIMIT = "TimeLimit";
constexpr auto DOT_PAR_OPTTOLERANCE = "OptTolerance";

constexpr auto DOT_PAR_RECODE = "Recode";

// My Network Simplex
#include "DOT_NetSimplex.h"

struct coprimes_t {
public:
  coprimes_t(int _v, int _w, double _c) : v(_v), w(_w), c_vw(_c) {}
  int v;
  int w;
  double c_vw;
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

typedef std::unordered_map<int_pair, double, pair_hash> int_pair_dict;
typedef std::unordered_map<int_pair, int, pair_hash> intpair2int;
typedef std::unordered_set<int_pair, pair_hash> int_pair_set;

namespace DOT {

class Histogram2D {
public:
  // Standard c'tor
  Histogram2D() {}

  // Second c'tor
  Histogram2D(int n, int *X, int *Y, int64_t *W) {
    for (int i = 0; i < n; i++)
      update(X[i], Y[i], W[i]);
  }

  // Third c'tor
  Histogram2D(int _n) : n(_n) { data.resize(_n * _n, 0); }

  // parse data from file
  void parse(const std::string &filename, const char sep = ',') {
    std::ifstream in_file(filename);

    if (!in_file) {
      fprintf(stdout, "FATAL ERROR: Cannot open file %s.\n", filename.c_str());
      exit(EXIT_FAILURE);
    }

    // Read first line
    auto read_row = [&](size_t i) {
      size_t j = 0;
      std::string line;
      std::getline(in_file, line);
      std::stringstream lineStream(line);
      std::string cell;

      while (std::getline(lineStream, cell, sep)) {
        data.push_back(stoll(cell));
        ++j;
      }

      return j;
    };

    // Read first row, and return row length
    n = read_row(0);

    for (size_t i = 1; i < n; ++i)
      read_row(i);

    // Use as few memory as possible
    data.shrink_to_fit();

    in_file.close();
  }

  // Total Weigth
  int64_t balance() {
    int64_t t = 0;
    for (const auto &k : data)
      t += k;
    return t;
  }

  // Support for loops
  std::vector<int64_t>::iterator begin() { return data.begin(); }
  std::vector<int64_t>::const_iterator begin() const { return data.begin(); }

  std::vector<int64_t>::iterator end() { return data.end(); }
  std::vector<int64_t>::const_iterator end() const { return data.end(); }

  // Add a new point
  void add(int _x, int _y, int64_t _w) { data[_x * n + _y] = _w; }

  // Add or update a new point
  void update(int _x, int _y, int64_t _w) { data[_x * n + _y] = _w; }

  // Getters
  size_t size() const { return data.size(); }

  // Get an elelment. Index data as a matrix
  int64_t get(size_t x, size_t y) const { return data[x * n + y]; };

  // Get dimension of histogram
  int getN() const { return n; }

private:
  int n; // Histogram of size n*n
  std::vector<int64_t>
      data; // Histogram data a single array (contiguos in memory)
};          // namespace DOT

// Solver class, which wrapper the Network Simplex algorithm
class Solver {
public:
  // Standard c'tor
  Solver()
      : _runtime(0.0), _n_log(0), L(-1), verbosity(DOT_VAL_INFO), recode(""),
        opt_tolerance(1e-06), timelimit(std::numeric_limits<double>::max()) {}

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

  // Compute Kantorovich-Wasserstein distance between two measures
  double solve(const Histogram2D &A, const Histogram2D &B) {
    auto start_t = std::chrono::steady_clock::now();
    double _all_p = 0.0;

    size_t n = A.getN();

    // Compute distances
    std::set<int64_t> tauset;
    for (size_t v = 0; v < n; ++v)
      for (size_t w = 0; w < n; ++w)
        tauset.insert(pow(v, 2) + pow(w, 2));

    vector<int64_t> tau;
    for (auto v : tauset)
      tau.push_back(v);

    fprintf(stdout, "distances: %d\n", tau.size());

    size_t idxL = 0;
    init_coprimes(tau[idxL++]);

    typedef int64_t FlowType;
    typedef int64_t CostType;

    // Build the graph for min cost flow
    NetSimplex<FlowType, CostType> simplex('E', static_cast<int>(2 * n * n), 0);

    // Set the parameters
    simplex.setTimelimit(timelimit);
    simplex.setVerbosity(verbosity);
    simplex.setOptTolerance(opt_tolerance);

    auto ID = [&n](int x, int y) { return x * n + y; };

    // add first d source nodes
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        simplex.addNode(ID(i, j), A.get(i, j));

    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        simplex.addNode(n * n + ID(i, j), -B.get(i, j));

    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j) {
        for (const auto &p : coprimes) {
          int v = p.v;
          int w = p.w;
          if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
            simplex.addArc(ID(i, j), n * n + ID(i + v, j + w), p.c_vw);
          }
        }
      }

    int it = 0;
    int n_cuts = 0;
    CostType fobj = 0;

    vector<double> pi(n, 0);

    Vars vars(n);
    for (size_t i = 0; i < n; ++i)
      vars[i].a = int(i);

    Vars vnew;
    vnew.reserve(n);

    // Init the simplex
    simplex.run();
    _iterations = simplex.iterations();

    // Start separation
    while (true) {
      _status = simplex.reRun();

      if (_status == ProblemType::TIMELIMIT)
        break;

      // Check feasibility: if feasible stop
      auto dummyFlow = simplex.dummyFlow();
      fprintf(stdout, "Flow: %ld, idx: %d\n", dummyFlow, idxL);

      if (abs(dummyFlow) < 1e-06 || idxL >= tau.size())
        break;

      // Update cost for dummy arcs (!)
      // MA SERVE? : simplex.updateDummyCost(tau[idxL+1]);

      // Add arcs
      init_coprimes(tau[idxL++]);

      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) {
          for (const auto &p : coprimes) {
            int v = p.v;
            int w = p.w;
            if (i + v >= 0 && i + v < n && j + w >= 0 && j + w < n) {
              simplex.addArc(ID(i, j), n * n + ID(i + v, j + w), p.c_vw);
            }
          }
        }

      ++it;
    }

    auto dummyFlow = simplex.dummyFlow();
    fprintf(stdout, "Flow: %ld, idx: %d\n", dummyFlow, idxL);

    _runtime = simplex.runtime();
    _iterations = simplex.iterations();
    _num_arcs = simplex.num_arcs();
    _num_nodes = simplex.num_nodes();

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    fobj = simplex.totalCost();

    if (_n_log > 0)
      PRINT("it: %d, fobj: %f, all: %f, simplex: %f, all_p: %f\n", it, fobj,
            _all, _runtime, _all_p);

    return fobj;
  }

  void init_coprimes(int64_t L) {
    int64_t L2 = L * L;
    coprimes.clear();
    for (int v = -L; v <= L; ++v)
      for (int w = -L; w <= L; ++w) {
        if (!(v == 0 && w == 0) && pow(v, 2) + pow(w, 2) == L2)
          coprimes.emplace_back(v, w, 1); // L2);
      }
    coprimes.shrink_to_fit();
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

  // List of pair of coprimes number between (-L, L)
  std::vector<coprimes_t> coprimes;

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

}; // namespace KWD

} // namespace DOT
