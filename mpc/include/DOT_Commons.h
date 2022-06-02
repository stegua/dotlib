/**
 * @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
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

#ifdef MY_RCPP
#include <R_ext/Print.h>
#define PRINT Rprintf
#else
#define PRINT printf
#endif  // MY_RCPP

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

void tolower(std::string& data) {
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

// NAME SPACE
namespace DOT {
enum class ProblemType {
  INFEASIBLE = 0,
  OPTIMAL = 1,
  UNBOUNDED = 2,
  TIMELIMIT = 3
};

enum class PivotRule { BLOCK_SEARCH = 0 };

class Var {
 public:
  int a;  // First point
  int b;  // Second point
  int c;  // Distance

  Var() : a(0), b(0), c(-1) {}

  Var(int _a, int _b, int _c) : a(_a), b(_b), c(_c) {}
};

typedef std::vector<Var> Vars;

class CapVar {
 public:
  int a;     // First point
  int b;     // Second point
  double c;  // Cost
  double d;  // Capacity

  CapVar() : a(0), b(0), c(-1), d(-1) {}

  CapVar(int _a, int _b, double _c, double _d) : a(_a), b(_b), c(_c), d(_d) {}
};

typedef std::vector<CapVar> CapVars;

// Read a text file and store it in a unique string
std::string readTextFile(const std::string& filename, size_t& len) {
  // Adapted from: https://github.com/simdjson/simdjson
  std::FILE* fp = std::fopen(filename.data(), "rb");

  if (fp == nullptr) return "";

  // Get the file size
  if (std::fseek(fp, 0, SEEK_END) < 0) {
    std::fclose(fp);
    return "";
  }

  long llen = std::ftell(fp);
  if ((llen < 0) || (llen == LONG_MAX)) {
    std::fclose(fp);
    return "";
  }

  // Allocate the string
  len = static_cast<size_t>(llen);
  std::string s;
  s.reserve(len);
  if (s.data() == nullptr) {
    std::fclose(fp);
    return "";
  }

  // Read the string
  std::rewind(fp);
  size_t bytes_read = std::fread((char*)(s.data()), 1, len, fp);
  if (std::fclose(fp) != 0 || bytes_read != len) return "";

  return s;
}

}  // namespace DOT