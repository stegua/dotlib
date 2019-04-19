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

// Compute Kantorovich-Wasserstein distance between two measures defined in R2 (bipartite version)
void BipartiteNetworkSimplex(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg);

void BipartiteCostScaling(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg);

void BipartiteCycleCanceling(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg);

void BipartiteEMD(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg);

} // End namespace
