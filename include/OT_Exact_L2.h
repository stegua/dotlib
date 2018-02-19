/**
* @fileoverview Copyright (c) 2017-2018, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "OT_BasicTypes.h"

/**
* @brief Solve standard Trasportation Problem on the complete bipartite graph
*        using Wassertein distance W_1^2 (order 1, with cost=L_2)
*/
real_t solve_exact_W_1_2(const histogram_t& h1, const histogram_t& h2);

real_t solve_exact_W_1_1(const histogram_t& h1, const histogram_t& h2);

real_t solve_exact_W_1_8(const histogram_t& h1, const histogram_t& h2);