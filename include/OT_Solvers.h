/**
* @fileoverview Copyright (c) 2017, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "OT_BasicTypes.h"

namespace DOTLib {
/**
* @brief Compute Earth Moving Distance (EMD) with L1 norm (with no power q)
*/
real_t solve_L1_1(const histogram_t& h1, const histogram_t& h2);

/**
* @brief Compute Earth Moving Distance (EMD) with L_infinity norm (with no power q)
*/
real_t solve_Linf_1(const histogram_t& h1, const histogram_t& h2);

};
