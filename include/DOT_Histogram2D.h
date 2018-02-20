/**
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "DOT_BasicTypes.h"

#include <fstream>
#include <sstream>

namespace DOT {

/**
* @brief Two dimensional histogram
*/
class Histogram2D {
 public:
   // Std c'tor
   Histogram2D(const std::string& filename) {
      readFromFile(filename);
   }

   // Rule of five?
   // ...

   // Parse from file
   void readFromFile(const std::string& filename) {
      std::ifstream in_file(filename);

      if (!in_file) {
         fprintf(stdout, "FATAL ERROR: Cannot open file %s.\n", filename.c_str());
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
            data.push_back(stoll(cell));
            ++j;
         }

         return j;
      };

      // Read first row, and return row length
      n = read_row(0);

      for (size_t i = 1; i < n; ++i)
         read_row(i);

      // Use as few memory as possible
      data.shrink_to_fit();

      in_file.close();
   }

   // Get an elelment. Index data as a matrix
   int64_t get(size_t x, size_t y) const {
      return data[x * n + y];
   };

   // Get dimension of histogram
   size_t getN() const {
      return n;
   }

   // Dump file to stdout
   void dump() const {
      for (size_t i = 0; i < n; ++i) {
         for (size_t j = 0; j < n; ++j)
            fprintf(stdout, "%ld ", get(i, j));
         fprintf(stdout, "\n");
      }
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
};

} // End namespace
