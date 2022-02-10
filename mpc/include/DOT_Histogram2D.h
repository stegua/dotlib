/**
 * @fileoverview Copyright (c) 2019-2022, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <DOT_Commons.h>

namespace DOT {

class Histogram2D {
 public:
  // Standard c'tor
  Histogram2D() {}

  // copy c'tor
  Histogram2D(const Histogram2D& o)
      : n(o.n), data(o.data.begin(), o.data.end()) {}

  // Third c'tor
  Histogram2D(int _n) : n(_n) { data.resize(_n * _n, 0); }

  // parse data from file
  void parse(const std::string& filename, const char sep = ',') {
    std::ifstream in_file(filename);

    if (!in_file) {
      fprintf(stdout, "FATAL ERROR: Cannot open file %s.\n", filename.c_str());
      exit(EXIT_FAILURE);
    }

    // Read first line
    auto read_row = [&](size_t i) {
      int j = 0;
      std::string line;
      std::getline(in_file, line);
      std::stringstream lineStream(line);
      std::string cell;

      while (std::getline(lineStream, cell, sep)) {
        data.push_back(stoi(cell));
        ++j;
      }

      return j;
    };

    // Read first row, and return row length
    n = read_row(0);

    for (size_t i = 1; i < n; ++i) read_row(i);

    // Use as few memory as possible
    data.shrink_to_fit();

    in_file.close();
  }

  // Total Weigth
  int64_t balance() const {
    int64_t t = 0;
    for (const auto& k : data) t += k;
    return t;
  }

  // Support for loops
  std::vector<int>::iterator begin() { return data.begin(); }
  std::vector<int>::const_iterator begin() const { return data.begin(); }

  std::vector<int>::iterator end() { return data.end(); }
  std::vector<int>::const_iterator end() const { return data.end(); }

  // Add a new point
  void add(size_t _x, size_t _y, int _w) { data[_x * n + _y] = _w; }

  // Add or update a new point
  void update(size_t _x, size_t _y, int _w) { data[_x * n + _y] = _w; }

  // Getters
  size_t size() const { return data.size(); }

  // Get an elelment. Index data as a matrix
  int get(size_t x, size_t y) const { return data[x * n + y]; };

  void subtract(size_t x, size_t y, int value) { data[x * n + y] -= value; };

  void set(size_t x, size_t y, int value) { data[x * n + y] = value; };

  // Get dimension of histogram
  int getN() const { return n; }

 private:
  // Histogram of size n*n
  int n;
  // Histogram data a single array (contiguos in memory)
  std::vector<int> data;
};

}  // namespace DOT
