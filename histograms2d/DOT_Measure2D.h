/**
 * @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#pragma once

#include <vector>
using std::vector;

#include <cmath>


namespace DOT {
class Measure2D {
 public:
  // Standard c'tor
  Measure2D() {}

  // Setter
  void reserve(size_t n) {
    Xs.reserve(n);
    Ys.reserve(n);
    Ws.reserve(n);
  }

  // Add a new point
  void add(size_t _x, size_t _y, double _w) {
    Xs.emplace_back(_x);
    Ys.emplace_back(_y);
    Ws.emplace_back(_w);
  }

  // Use as few memory as possible
  void shrink_to_fit() {
    Xs.shrink_to_fit();
    Ys.shrink_to_fit();
    Ws.shrink_to_fit();
  }

  // Getters
  size_t size() const { return Ws.size(); }

  inline double getWeight(size_t i) const { return Ws[i]; }

  inline double getX(size_t i) const { return Xs[i]; }
  inline double getY(size_t i) const { return Ys[i]; }

 private:
  vector<size_t> Xs;
  vector<size_t> Ys;
  vector<double> Ws;
};


class Plan {
public:
  // Setter
  void reserve(size_t n, size_t m) {
    head.reserve(n + m - 1);
    tail.reserve(n + m - 1);
    flow.reserve(n + m - 1);
  }

  // Add a new point
  void add(size_t _i, size_t _j, double _w) {
    tail.emplace_back(_i);
    head.emplace_back(_j);
    flow.emplace_back(_w);
  }

  // Use as few memory as possible
  void shrink_to_fit() {
    tail.shrink_to_fit();
    head.shrink_to_fit();
    flow.shrink_to_fit();
  }

  // Getters
  size_t size() const { return flow.size(); }

private:
    // Amount of flow in arc (tail, head)
    vector<size_t> tail;
    vector<size_t> head;
    vector<double> flow;  
};


// Solve Optimal Transport between two 1D measures
//std::pair<double, Plan> Wasserstein2_1D(const Measure1D& mu, const Measure1D& nu);

} // end namespace DOT
