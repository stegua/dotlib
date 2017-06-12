/**
* @fileoverview Copyright (c) 2017, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

// FROM STDLIB
#include <cassert>

// For PI constant definition as M_PI
#define _USE_MATH_DEFINES
#include <math.h>

// My Types for OT lib
typedef double   real_t;

#include <vector>
typedef std::vector< real_t >  histogram_t;
typedef std::vector< real_t >  row_t;
typedef std::vector< row_t  >  matrix_t;

#include <unordered_set>
typedef std::unordered_set<int> int_set;

struct pair_hash {
   inline std::size_t operator()(const std::pair<int, int> & v) const {
      return v.first * 31 + v.second;
   }
};
typedef std::unordered_set<std::pair<int, int>, pair_hash> pair_set;


// LEMON GRAPH LIBRARY
#include <lemon/smart_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <lemon/cycle_canceling.h>

#include <lemon/preflow.h>

using namespace lemon;
typedef lemon::ListDigraph Graph;
typedef int64_t LimitValueType;


// EXTERNAL LIBRARIES:

// SpdLog (https://github.com/gabime/spdlog)
#include "spdlog/spdlog.h"
namespace spd = spdlog;
