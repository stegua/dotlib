/**
* @fileoverview Copyright (c) 2017, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "OT_BasicTypes.h"

/**
* @brief Compute Earth Moving Distance (EMD) with L1 norm (with no power)
*        between pair of images spacing razs at "constant" angle distance, with given degree ration
*        Example: with degree=45, it puts 8 rays for vertex (at most)
*/
real_t solve_L1_1(const histogram_t& h1, const histogram_t& h2);
