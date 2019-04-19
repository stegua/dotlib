/*
*  Main authors:
*     @author Stefano Gualandi <stefano.gualandi@gmail.com>
*
*     @fileoverview Copyright (c) 2017-19, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
*  Last update: April, 2019
*/

#pragma once


#include "DOT_MeasureR2.hpp"

namespace DOT {

// Compute EXACT Kantorovich-Wasserstein Distance (KWD) of order 1,
// using the BGV method with ground distance L1 (Manhattan)
void BGV_L1(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg);

// Compute EXACT Kantorovich-Wasserstein Distance (KWD) of order 1,
// using the BGV method with ground distance L_inf
void BGV_Linf(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg);

// Compute EXACT Kantorovich-Wasserstein Distance (KWD) of order 1,
// using the BGV method with ground distance L2
void BGV_L2(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg,
            const std::vector<std::pair<int, int>>& coprimes);



// Compute Kantorovich-Wasserstein Distance (KWD) of order 1, with ground distance L1 (Manhattan)
// using the original code of Link&Okada (see externs/ source files)
void EMD_L1(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg);

} // End namespace
