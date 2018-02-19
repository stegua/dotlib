/**
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

// Best debbuger ever
#define KABOOM(CONDITION, MSG)           \
	if(CONDITION) {                      \
		fprintf(stdout, "KA KA KA BOOM: %s!!!\n", MSG);  \
		exit(EXIT_FAILURE);			     \
    }


// FROM STDLIB
#include <cassert>
#include <string>
#include <vector>

// For PI constant definition as M_PI
#define _USE_MATH_DEFINES
#include <math.h>


// LEMON GRAPH LIBRARY
#include <lemon/smart_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <lemon/cost_scaling.h>
#include <lemon/cycle_canceling.h>

#include <lemon/preflow.h>

using namespace lemon;
typedef lemon::ListDigraph LemonGraph;
typedef int64_t LimitValueType;
typedef NetworkSimplex<LemonGraph, LimitValueType, LimitValueType> LemonSimplex;
