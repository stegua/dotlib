/**
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "yocta_logger.hh"

#include "DOT_BasicTypes.h"

namespace DOT {

/**
* @brief Two dimensional histogram
*/
class Histogram {
 public:
   // Std c'tor
   Histogram(const std::string& filename, size_t _n, size_t _d) : n(_n), d(_d) {
      readFromCSVFile(filename);
      build();
   }

   // Rule of five?
   // ...

   // Parse from file
   void readFromCSVFile(const std::string& filename) {
      std::ifstream in_file(filename);

      if (!in_file) {
         logger.error("FATAL ERROR: Cannot open file %s", filename.c_str());
         exit(EXIT_FAILURE);
      }

      std::string csv_file = yocta::readTextFile(filename);
      std::istringstream input;
      input.str(csv_file);
      std::istringstream in_row;
      std::vector<float> tmp(d, 0);
      bool header = true;
      for (std::string row; std::getline(input, row); ) {
         if (header) {
            header = false;
         } else {
            tmp.clear();

            in_row.clear();
            in_row.str(row);
            int i = 0;
            // Specific for FCS files (TODO: refactor)
            for (std::string col; std::getline(in_row, col, ','); i++) {
               if (i == 1) {
                  tmp.push_back(std::stof(col.c_str()));
               } else {
                  if (i > 1 && i <= d) {
                     float t = std::stof(col.c_str());
                     tmp.push_back(std::pow(10.0, t/ 256));
                  }
               }
            }

            raw_data.emplace_back(tmp);
         }
      }
   }

   // Build a d-dimensional histogram
   void build() {
      // 2-dimensional
      if (d == 2) {
         data.resize(n*n, 0);
         for (const auto& v : raw_data) {
            assert(v.size() == 2);
            size_t x = std::max<float>(0.0, std::min<float>(n-1, floor(v[0] / 1024 * n)));
            size_t y = std::max<float>(0.0, std::min<float>(n-1, floor(v[1] / 1024 * n)));

            data[x*n + y]++;
         }
      }

      // 3-dimensional
      if (d == 3) {
         data.resize(n*n*n, 0);
         for (const auto& v : raw_data) {
            assert(v.size() == 3);
            size_t x = std::max<float>(0.0, std::min<float>(n - 1, floor(v[0] / 1024 * n)));
            size_t y = std::max<float>(0.0, std::min<float>(n - 1, floor(v[1] / 1024 * n)));
            size_t z = std::max<float>(0.0, std::min<float>(n - 1, floor(v[2] / 1024 * n)));

            data[x*n*n + y*n + z]++;
         }
      }

      // 4-dimensional
      if (d == 4) {
         data.resize(n*n*n*n, 0);
         for (const auto& v : raw_data) {
            assert(v.size() == 4);
            size_t x = std::max<float>(0.0, std::min<float>(n - 1, floor(v[0] / 1024 * n)));
            size_t y = std::max<float>(0.0, std::min<float>(n - 1, floor(v[1] / 1024 * n)));
            size_t z = std::max<float>(0.0, std::min<float>(n - 1, floor(v[2] / 1024 * n)));
            size_t w = std::max<float>(0.0, std::min<float>(n - 1, floor(v[3] / 1024 * n)));

            data[x*n*n*n + y * n*n + z*n + w]++;
         }
      }
   }

   // Get an elelment. Index data as a matrix
   int64_t get(size_t x, size_t y) const {
      return data[x * n + y];
   };
   int64_t get(size_t x, size_t y, size_t z) const {
      return data[x*n*n + y * n + z];
   };
   int64_t get(size_t x, size_t y, size_t z, size_t w) const {
      return data[x*n*n*n + y * n*n + z * n + w];
   };

   // Get dimension of histogram
   size_t getN() const {
      return n;
   }

   // Get dimension of histogram
   size_t getD() const {
      return d;
   }

   // Dump file to stdout
   void dump() const {
   }

   // Sum over all bins
   int64_t computeTotalWeight() const {
      int64_t t = 0;
      for (int64_t v : data)
         t += v;
      return t;
   }

 private:
   size_t n;   // Histogram of size n*n
   std::vector<int64_t>  data;   // Histogram data a single array (contiguos in memory)

   size_t d;  // dimension of the histogram
   std::vector<std::vector<float>> raw_data; // raw data, one row for each single point
};

} // End namespace
