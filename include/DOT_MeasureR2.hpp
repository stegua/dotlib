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

#include <array>
using std::array;

#include <vector>
using std::vector;

#include <fstream>
#include <sstream>


namespace DOT {
// Integer tolerance for double distances
const double INTEGER_TOL = 1000000;

// Distance function between a pair of point in R^k
#define DISTANCE_R2(x, y)  std::round(sqrt((x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]))*INTEGER_TOL)
#define DISTANCE_R1(x, y)  fabs(x[0] - y[0]) + fabs(x[1] - y[1])
#define DISTANCE_Rinf(x, y)  std::max<>(fabs(x[0] - y[0]), fabs(x[1] - y[1]))

//#define DISTANCE_R2(x, y)  (x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1])


// Discrete Optimal Transport in points of R^2
typedef std::array<int, 2> PointR2;

// Container for general discrete measure
class MeasureR2 {
 public:
   MeasureR2() {}

   MeasureR2(const std::string& filename) {
      readFromFile(filename);
   }

   // setter
   void reserve(size_t n) {
      Ws.reserve(n);
      Ps.reserve(n);
   }

   void add(float _w, PointR2 _p) {
      Ws.emplace_back(_w);
      Ps.emplace_back(_p);
   }

   // getters
   size_t size() const {
      return Ws.size();
   }

   vector<double> getWs() const {
      return Ws;
   }

   float getW(size_t i) const {
      return Ws[i];
   }

   vector<PointR2> getPs() const {
      return Ps;
   }

   PointR2 getP(size_t i) const {
      return Ps[i];
   }

   // Parse from file
   void readFromFile(const std::string& filename) {
      std::ifstream in_file(filename);

      if (!in_file) {
         fprintf(stdout, "FATAL ERROR: Cannot open file %s", filename.c_str());
         exit(EXIT_FAILURE);
      }

      // Read first line
      auto read_row = [&](size_t i) {
         size_t j = 0;
         std::string         line;
         std::getline(in_file, line);
         std::stringstream   lineStream(line);
         std::string         cell;

         while (std::getline(lineStream, cell, ',')) {
            Ws.push_back(stof(cell));
            Ps.emplace_back(PointR2({ int(i), int(j) }));
            ++j;
         }

         return j;
      };

      // Read first row, and return row length
      auto n = read_row(0);

      Ws.reserve(n*n);
      Ps.reserve(n*n);

      for (size_t i = 1; i < n; ++i)
         read_row(i);

      in_file.close();

      // Use as few memory as possible
      Ws.shrink_to_fit();
      Ps.shrink_to_fit();
   }

 private:
   vector<double>    Ws;
   vector<PointR2>  Ps;
};

// Create a random measure in R2
MeasureR2 createRandom0N(size_t n, int seed);

// Compute co-prime numbers over a lattice n*n with max size L
// List of pair of coprimes number between (-L, L)
std::vector<std::pair<int, int>> buildCoprimes(int L);

} // End namespace
