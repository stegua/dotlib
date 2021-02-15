/**
 * @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */


#pragma once

#include <chrono>
#include <exception>
#include <limits>
#include <vector>
#include <string>

#ifdef MY_RCPP
#include <R_ext/Print.h>
#define PRINT Rprintf
#else
#define PRINT printf
#endif // MY_RCPP

namespace DOT {

	const double FEASIBILITY_TOL = 1e-09;
	const double PRIC_TOL = 1e-09;

	enum class ProblemType {
		INFEASIBLE = 0,
		OPTIMAL = 1,
		UNBOUNDED = 2,
		TIMELIMIT = 3
	};

	enum class PivotRule { BLOCK_SEARCH = 0 };

	template <typename V = int, typename C = V> class GVar {
	public:
		V a; // First point
		V b; // Second point
		C c; // Distance

		GVar() : a(0), b(0), c(-1) {}

		GVar(V _a, V _b, C _c) : a(_a), b(_b), c(_c) {}
	};

	typedef GVar<int, double> Var;
	typedef std::vector<Var> Vars;

} // END namespace