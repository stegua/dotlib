/**
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "DOT_BasicTypes.h"
#include "DOT_Histogram2D.h"
#include "DOT_NetSimplex.h"

namespace DOT {

	/**
	* @brief Compute the Wasserstein distance of order 1 with ground distance L1
	*/
	double solve_network_W1L2(
		const Histogram2D& h1,
		const Histogram2D& h2,
		std::vector<std::pair<int, int>> coprimes)
	{

	}  // end namespace DOT
