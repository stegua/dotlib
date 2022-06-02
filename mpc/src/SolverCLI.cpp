/*
 * @fileoverview Copyright (c) 2019-2021, Stefano Gualandi,
 *               via Ferrata, 5, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#include <algorithm>
#include <cassert>
#include <random>

#include "DOT_CoinLemon.h"
#include "DOT_Cplex.h"
#include "DOT_EapiSimplex.h"
#include "DOT_Lemon.h"
#include "DOT_NetSimplex.h"
#include "DOT_Solver.h"

// Zeta block for coordinates vector (CUDA
#define BLOCKSIZE 1024
#define BLOCKNUM 8

int seed = 13;

std::random_device
    rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(seed);  // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> Uniform01(0, 1);

int testEapi1() {
  DOT::EapiSimplex solver(4, 4);

  solver.addNode(0, 1);
  solver.addNode(1, 1);
  solver.addNode(2, -1);
  solver.addNode(3, -1);

  solver.addArc(0, 2, 1);
  solver.addArc(0, 3, 4);
  solver.addArc(1, 2, 0);
  solver.addArc(1, 3, 1);

  int r = solver.solve();
  return r;
}

int testEapi2() {
  DOT::EapiSimplex solver(6, 9);

  //  0   1   2
  //    \   \   \
  //      3   4   5
  solver.addNode(0, 2);
  solver.addNode(1, 1);
  solver.addNode(2, 0);

  solver.addNode(3, -0);
  solver.addNode(4, -1);
  solver.addNode(5, -2);

  solver.addArc(0, 3, 1);
  solver.addArc(0, 4, 4);
  solver.addArc(0, 5, 9);

  solver.addArc(1, 3, 0);
  solver.addArc(1, 4, 1);
  solver.addArc(1, 5, 4);

  solver.addArc(2, 3, 1);
  solver.addArc(2, 4, 0);
  solver.addArc(2, 5, 1);

  int r = solver.solve();
  return r;
}

int testEapi3() {
  DOT::EapiSimplex solver(4, 8);

  //  0 - 1
  //  |   |
  //  2 - 3

  solver.addNode(0, 1);
  solver.addNode(1, 0);
  solver.addNode(2, 0);
  solver.addNode(3, -1);

  solver.addArc(0, 1, 1);
  solver.addArc(1, 0, 1);
  solver.addArc(0, 2, 1);
  solver.addArc(2, 0, 1);

  solver.addArc(3, 1, 1);
  solver.addArc(1, 3, 1);
  solver.addArc(3, 2, 1);
  solver.addArc(2, 3, 1);

  int r = solver.solve();

  return r;
}

int testEapi4() {
  DOT::EapiSimplex solver(6, 7);

  //  0 - 1 - 2
  //  |   |   |
  //  3 - 4 - 5

  solver.addNode(0, 1);
  solver.addNode(1, 0);
  solver.addNode(2, 0);

  solver.addNode(3, 0);
  solver.addNode(4, 0);
  solver.addNode(5, -1);

  solver.addArc(0, 1, 1);  // 7
  solver.addArc(0, 3, 1);  // 8
  solver.addArc(1, 2, 1);  // 9
  solver.addArc(1, 4, 1);  // 10
  solver.addArc(2, 5, 1);  // 11
  solver.addArc(3, 4, 1);  // 12
  solver.addArc(4, 5, 1);  // 13

  int r = solver.solve();
  return r;
}

int testEapi5() {
  DOT::EapiSimplex solver(9, 12);

  //  0 - 1 - 2
  //  |   |   |
  //  3 - 4 - 5
  //  |   |   |
  //  6 - 7 - 8

  solver.addNode(0, 2);
  solver.addNode(1, 1);
  solver.addNode(2, 0);

  solver.addNode(3, 0);
  solver.addNode(4, 0);
  solver.addNode(5, 0);

  solver.addNode(6, 0);
  solver.addNode(7, -1);
  solver.addNode(8, -2);

  solver.addArc(0, 1, 1);  // 9
  solver.addArc(0, 3, 1);  // 10
  solver.addArc(1, 2, 1);  // 11
  solver.addArc(1, 4, 1);  // 12
  solver.addArc(2, 5, 1);  // 13
  solver.addArc(3, 4, 1);  // 14
  solver.addArc(3, 6, 1);  // 15
  solver.addArc(4, 5, 1);  // 16
  solver.addArc(4, 7, 1);  // 17
  solver.addArc(5, 8, 1);  // 18
  solver.addArc(6, 7, 1);  // 19
  solver.addArc(7, 8, 1);  // 20

  int r = solver.solve();
  return r;
}

int testEapi6() {
  DOT::EapiSimplex solver(9, 12);

  //  0 - 1 - 2
  //  |   |   |
  //  3 - 4 - 5
  //  |   |   |
  //  6 - 7 - 8

  solver.addNode(0, 2);
  solver.addNode(1, -1);
  solver.addNode(2, 0);

  solver.addNode(3, 0);
  solver.addNode(4, 0);
  solver.addNode(5, 0);

  solver.addNode(6, 0);
  solver.addNode(7, +1);
  solver.addNode(8, -2);

  solver.addArc(0, 1, 1);  // 9
  solver.addArc(0, 3, 1);  // 10
  solver.addArc(1, 2, 1);  // 11
  solver.addArc(1, 4, 1);  // 12
  solver.addArc(2, 5, 1);  // 13
  solver.addArc(3, 4, 1);  // 14
  solver.addArc(3, 6, 1);  // 15
  solver.addArc(4, 5, 1);  // 16
  solver.addArc(4, 7, 1);  // 17
  solver.addArc(5, 8, 1);  // 18
  solver.addArc(6, 7, 1);  // 19
  solver.addArc(7, 8, 1);  // 20

  int r = solver.solve();
  return r;
}

int testEapi7() {
  DOT::EapiSimplex solver(9, 24);

  //  0 - 1 - 2
  //  |   |   |
  //  3 - 4 - 5
  //  |   |   |
  //  6 - 7 - 8

  solver.addNode(0, 3);
  solver.addNode(1, 0);
  solver.addNode(2, -2);

  solver.addNode(3, 0);
  solver.addNode(4, 0);
  solver.addNode(5, 0);

  solver.addNode(6, -1);
  solver.addNode(7, -1);
  solver.addNode(8, +1);

  solver.addArc(0, 1, 1);  // 9
  solver.addArc(0, 3, 1);  // 10
  solver.addArc(1, 2, 1);  // 11
  solver.addArc(1, 4, 1);  // 12
  solver.addArc(2, 5, 1);  // 13
  solver.addArc(3, 4, 1);  // 14
  solver.addArc(3, 6, 1);  // 15
  solver.addArc(4, 5, 1);  // 16
  solver.addArc(4, 7, 1);  // 17
  solver.addArc(5, 8, 1);  // 18
  solver.addArc(6, 7, 1);  // 19
  solver.addArc(7, 8, 1);  // 20

  solver.addArc(1, 0, 1);  // 21
  solver.addArc(3, 0, 1);  // 22
  solver.addArc(2, 1, 1);  // 23
  solver.addArc(4, 1, 1);  // 24
  solver.addArc(5, 2, 1);  // 25
  solver.addArc(4, 3, 1);  // 26
  solver.addArc(6, 3, 1);  // 27
  solver.addArc(5, 4, 1);  // 28
  solver.addArc(7, 4, 1);  // 29
  solver.addArc(8, 5, 1);  // 30
  solver.addArc(7, 6, 1);  // 31
  solver.addArc(8, 7, 1);  // 32

  int r = solver.solve();
  return r;
}

int testEapi8() {
  for (int n = 5; n < 30; n++) {
    std::uniform_int_distribution<> Uniform1N(1, n);

    DOT::EapiSimplex solver(2 * n, n * n);

    DOT::NetSimplex lemon('F', 2 * n, n * n);
    lemon.setTimelimit(3600);
    lemon.setOptTolerance(0.001);

    vector<int> tmp;
    for (size_t i = 0; i < n; i++) {
      auto d = Uniform1N(gen);
      tmp.push_back(d);

      solver.addNode(i, d);
      lemon.addNode(i, d);
    }

    std::shuffle(tmp.begin(), tmp.end(), gen);
    for (size_t i = 0; i < n; i++) {
      solver.addNode(n + i, -tmp[i]);
      lemon.addNode(n + i, -tmp[i]);
    }

    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++) {
        solver.addArc(i, n + j, (i - j) * (i - j));
        lemon.addArc(i, n + j, (i - j) * (i - j));
      }

    int r = solver.solve();

    lemon.run();
    int s = (int)lemon.totalCost();
    fprintf(stdout, "CHECK %d ==> eapi: %d, lemon: %d\n", n, r, s);
    assert(s == r);
  }

  return 1;
}

int testEapi9() {
  for (int n = 5; n < 65; n++) {
    std::uniform_int_distribution<> Uniform1N(0, 255);

    DOT::EapiSimplex solver(2 * n * n, n * n * n * n);

    DOT::NetSimplex lemon('F', 2 * n * n, n * n * n * n);
    lemon.setTimelimit(3600);
    lemon.setOptTolerance(0.001);

    vector<int> tmp;
    for (size_t i = 0; i < n * n; i++) {
      auto d = Uniform1N(gen);
      tmp.push_back(d);

      solver.addNode(i, d);
      lemon.addNode(i, d);
    }

    std::shuffle(tmp.begin(), tmp.end(), gen);
    for (size_t i = 0; i < n * n; i++) {
      solver.addNode(n * n + i, -tmp[i]);
      lemon.addNode(n * n + i, -tmp[i]);
    }

    auto ID = [&n](int x, int y) { return x * n + y; };

    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++) {
        for (size_t v = 0; v < n; v++)
          for (size_t w = 0; w < n; w++) {
            solver.addArc(ID(i, j), n * n + ID(v, w),
                          (i - v) * (i - v) + (j - w) * (j - w));
            lemon.addArc(ID(i, j), n * n + ID(v, w),
                         (i - v) * (i - v) + (j - w) * (j - w));
          }
      }

    int r = solver.solve();
    int api_it = (int)solver.iterations();
    double api_time = solver.runtime();
    double api_pri = solver.pricing_time();

    lemon.run();
    int s = (int)lemon.totalCost();
    int lem_it = (int)lemon.iterations();
    double lem_time = lemon.runtime();
    fprintf(
        stdout,
        "CHECK %d %d ==> eapi: %d, it: %d, time: %.3f, pricing: %.4f | lemon: "
        "%d, it: %d, time: %.3f, pricing: %.4f\n",
        n, solver._block_size, r, api_it, api_time, api_pri, s, lem_it,
        lem_time, lemon._time_pricing);
    assert(s == r);
  }

  return 1;
}

int main(int argc, char *argv[]) {
  int n = 8;

  if (argc > 1) n = atoi(argv[1]);

  std::uniform_int_distribution<> Uniform1N(1, n);

  if (true) {
    DOT::SolverCoinLemon solver;
    size_t len = 0;

    auto start_t = std::chrono::steady_clock::now();

    // PA_dense_SS_tract. PA_dense_CD_county, PA_dense_SS_county,
    // PA_dense_CD_tract
    solver.parseDIMACS(
        "C:\\Users\\Gualandi\\Data\\Districting\\CA_dense3_SH_tract.net", len);
    // solver.parseDIMACS(
    //    "C:\\Users\\Gualandi\\Data\\Districting\\PA_dense_CD_county.net",
    //    len);

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    fprintf(stdout, "start solver: %f\n", _all);
    fflush(stdout);

    // solver.dump();
    solver.solve();
  }

  if (true) {
    DOT::SolverNSC solver;
    size_t len = 0;

    auto start_t = std::chrono::steady_clock::now();

    // PA_dense_SS_tract. PA_dense_CD_county, PA_dense_SS_county,
    // PA_dense_CD_tract
    solver.parseDIMACS(
        "C:\\Users\\Gualandi\\Data\\Districting\\CA_dense3_SH_tract.net", len);

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    fprintf(stdout, "start solver: %f\n", _all);
    fflush(stdout);

    // solver.dump();
    solver.colgen();
  }

  if (true) {
    DOT::SolverCplex solver;
    size_t len = 0;

    auto start_t = std::chrono::steady_clock::now();

    // PA_dense_SS_tract. PA_dense_CD_county, PA_dense_SS_county,
    // PA_dense_CD_tract
    parseDIMACS(
        "C:\\Users\\Gualandi\\Data\\Districting\\CA_dense3_SH_tract.net", len,
        solver);

    auto end_t = std::chrono::steady_clock::now();
    auto _all = double(std::chrono::duration_cast<std::chrono::milliseconds>(
                           end_t - start_t)
                           .count()) /
                1000;

    fprintf(stdout, "start solver: %f\n", _all);
    fflush(stdout);

    // solver.dump();
    solver.solve();
  }

  if (false) {
    assert(2 == testEapi1());
    assert(17 == testEapi2());
    assert(2 == testEapi3());
    assert(3 == testEapi4());
    assert(10 == testEapi5());
    assert(6 == testEapi6());

    // testEapi7();
    // testEapi8();
    //// testEapi9();
    // testEapi2();
  }

  if (false) {
    std::string SEP = "\\";

    std::string base = "C:\\Users\\Gualandi\\Data\\ColorTransfer";
    // std::vector<std::string> Fs = {"color-8192-seed-19.net"};
    // std::vector<std::string> Fs = {"color-4096.net"};

    std::vector<std::string> Fs = {
        "color-1024-seed-29.net",  "color-11264-seed-29.net",
        "color-12800-seed-29.net", "color-14336-seed-29.net",
        "color-2048-seed-29.net",  "color-3584-seed-29.net",
        "color-512-seed-29.net",   "color-6144-seed-29.net",
        "color-7680-seed-29.net",  "color-9216-seed-29.net",
        "color-10240-seed-29.net", "color-11776-seed-29.net",
        "color-13312-seed-29.net", "color-14848-seed-29.net",
        "color-2560-seed-29.net",  "color-4096-seed-29.net",
        "color-5120-seed-29.net",  "color-6656-seed-29.net",
        "color-8192-seed-29.net",  "color-9728-seed-29.net",
        "color-10752-seed-29.net", "color-12288-seed-29.net",
        "color-13824-seed-29.net", "color-1536-seed-29.net",
        "color-3072-seed-29.net",  "color-4608-seed-29.net",
        "color-5632-seed-29.net",  "color-7168-seed-29.net",
        "color-8704-seed-29.net",  "color-1024-seed-17.net",
        "color-11264-seed-17.net", "color-12800-seed-17.net",
        "color-14336-seed-17.net", "color-2048-seed-17.net",
        "color-3584-seed-17.net",  "color-512-seed-17.net",
        "color-6144-seed-17.net",  "color-7680-seed-17.net",
        "color-9216-seed-17.net",  "color-10240-seed-17.net",
        "color-11776-seed-17.net", "color-13312-seed-17.net",
        "color-14848-seed-17.net", "color-2560-seed-17.net",
        "color-4096-seed-17.net",  "color-5120-seed-17.net",
        "color-6656-seed-17.net",  "color-8192-seed-17.net",
        "color-9728-seed-17.net",  "color-10752-seed-17.net",
        "color-12288-seed-17.net", "color-13824-seed-17.net",
        "color-1536-seed-17.net",  "color-3072-seed-17.net",
        "color-4608-seed-17.net",  "color-5632-seed-17.net",
        "color-7168-seed-17.net",  "color-8704-seed-17.net",
        "color-1024-seed-19.net",  "color-11264-seed-19.net",
        "color-12800-seed-19.net", "color-14336-seed-19.net",
        "color-2048-seed-19.net",  "color-3584-seed-19.net",
        "color-512-seed-19.net",   "color-6144-seed-19.net",
        "color-7680-seed-19.net",  "color-9216-seed-19.net",
        "color-10240-seed-19.net", "color-11776-seed-19.net",
        "color-13312-seed-19.net", "color-14848-seed-19.net",
        "color-2560-seed-19.net",  "color-4096-seed-19.net",
        "color-5120-seed-19.net",  "color-6656-seed-19.net",
        "color-8192-seed-19.net",  "color-9728-seed-19.net",
        "color-10752-seed-19.net", "color-12288-seed-19.net",
        "color-13824-seed-19.net", "color-1536-seed-19.net",
        "color-3072-seed-19.net",  "color-4608-seed-19.net",
        "color-5632-seed-19.net",  "color-7168-seed-19.net",
        "color-8704-seed-19.net",  "color-1024-seed-23.net",
        "color-11264-seed-23.net", "color-12800-seed-23.net",
        "color-14336-seed-23.net", "color-2048-seed-23.net",
        "color-3584-seed-23.net",  "color-512-seed-23.net",
        "color-6144-seed-23.net",  "color-7680-seed-23.net",
        "color-9216-seed-23.net",  "color-10240-seed-23.net",
        "color-11776-seed-23.net", "color-13312-seed-23.net",
        "color-14848-seed-23.net", "color-2560-seed-23.net",
        "color-4096-seed-23.net",  "color-5120-seed-23.net",
        "color-6656-seed-23.net",  "color-8192-seed-23.net",
        "color-9728-seed-23.net",  "color-10752-seed-23.net",
        "color-12288-seed-23.net", "color-13824-seed-23.net",
        "color-1536-seed-23.net",  "color-3072-seed-23.net",
        "color-4608-seed-23.net",  "color-5632-seed-23.net",
        "color-7168-seed-23.net",  "color-8704-seed-23.net",
        "color-1024-seed-29.net",  "color-11264-seed-29.net",
        "color-12800-seed-29.net", "color-14336-seed-29.net",
        "color-2048-seed-29.net",  "color-3584-seed-29.net",
        "color-512-seed-29.net",   "color-6144-seed-29.net",
        "color-7680-seed-29.net",  "color-9216-seed-29.net",
        "color-10240-seed-29.net", "color-11776-seed-29.net",
        "color-13312-seed-29.net", "color-14848-seed-29.net",
        "color-2560-seed-29.net",  "color-4096-seed-29.net",
        "color-5120-seed-29.net",  "color-6656-seed-29.net",
        "color-8192-seed-29.net",  "color-9728-seed-29.net",
        "color-10752-seed-29.net", "color-12288-seed-29.net",
        "color-13824-seed-29.net", "color-1536-seed-29.net",
        "color-3072-seed-29.net",  "color-4608-seed-29.net",
        "color-5632-seed-29.net",  "color-7168-seed-29.net",
        "color-8704-seed-29.net",  "color-1024-seed-31.net",
        "color-11264-seed-31.net", "color-12800-seed-31.net",
        "color-14336-seed-31.net", "color-2048-seed-31.net",
        "color-3584-seed-31.net",  "color-512-seed-31.net",
        "color-6144-seed-31.net",  "color-7680-seed-31.net",
        "color-9216-seed-31.net",  "color-10240-seed-31.net",
        "color-11776-seed-31.net", "color-13312-seed-31.net",
        "color-14848-seed-31.net", "color-2560-seed-31.net",
        "color-4096-seed-31.net",  "color-5120-seed-31.net",
        "color-6656-seed-31.net",  "color-8192-seed-31.net",
        "color-9728-seed-31.net",  "color-10752-seed-31.net",
        "color-12288-seed-31.net", "color-13824-seed-31.net",
        "color-1536-seed-31.net",  "color-3072-seed-31.net",
        "color-4608-seed-31.net",  "color-5632-seed-31.net",
        "color-7168-seed-31.net",  "color-8704-seed-31.net",
        "color-1024-seed-37.net",  "color-11264-seed-37.net",
        "color-12800-seed-37.net", "color-14336-seed-37.net",
        "color-2048-seed-37.net",  "color-3584-seed-37.net",
        "color-512-seed-37.net",   "color-6144-seed-37.net",
        "color-7680-seed-37.net",  "color-9216-seed-37.net",
        "color-10240-seed-37.net", "color-11776-seed-37.net",
        "color-13312-seed-37.net", "color-14848-seed-37.net",
        "color-2560-seed-37.net",  "color-4096-seed-37.net",
        "color-5120-seed-37.net",  "color-6656-seed-37.net",
        "color-8192-seed-37.net",  "color-9728-seed-37.net",
        "color-10752-seed-37.net", "color-12288-seed-37.net",
        "color-13824-seed-37.net", "color-1536-seed-37.net",
        "color-3072-seed-37.net",  "color-4608-seed-37.net",
        "color-5632-seed-37.net",  "color-7168-seed-37.net",
        "color-8704-seed-37.net",  "color-1024-seed-41.net",
        "color-11264-seed-41.net", "color-12800-seed-41.net",
        "color-14336-seed-41.net", "color-2048-seed-41.net",
        "color-3584-seed-41.net",  "color-512-seed-41.net",
        "color-6144-seed-41.net",  "color-7680-seed-41.net",
        "color-9216-seed-41.net",  "color-10240-seed-41.net",
        "color-11776-seed-41.net", "color-13312-seed-41.net",
        "color-14848-seed-41.net", "color-2560-seed-41.net",
        "color-4096-seed-41.net",  "color-5120-seed-41.net",
        "color-6656-seed-41.net",  "color-8192-seed-41.net",
        "color-9728-seed-41.net",  "color-10752-seed-41.net",
        "color-12288-seed-41.net", "color-13824-seed-41.net",
        "color-1536-seed-41.net",  "color-3072-seed-41.net",
        "color-4608-seed-41.net",  "color-5632-seed-41.net",
        "color-7168-seed-41.net",  "color-8704-seed-41.net",
        "color-1024-seed-43.net",  "color-11264-seed-43.net",
        "color-12800-seed-43.net", "color-14336-seed-43.net",
        "color-2048-seed-43.net",  "color-3584-seed-43.net",
        "color-512-seed-43.net",   "color-6144-seed-43.net",
        "color-7680-seed-43.net",  "color-9216-seed-43.net",
        "color-10240-seed-43.net", "color-11776-seed-43.net",
        "color-13312-seed-43.net", "color-14848-seed-43.net",
        "color-2560-seed-43.net",  "color-4096-seed-43.net",
        "color-5120-seed-43.net",  "color-6656-seed-43.net",
        "color-8192-seed-43.net",  "color-9728-seed-43.net",
        "color-10752-seed-43.net", "color-12288-seed-43.net",
        "color-13824-seed-43.net", "color-1536-seed-43.net",
        "color-3072-seed-43.net",  "color-4608-seed-43.net",
        "color-5632-seed-43.net",  "color-7168-seed-43.net",
        "color-8704-seed-43.net",  "color-1024-seed-47.net",
        "color-11264-seed-47.net", "color-12800-seed-47.net",
        "color-14336-seed-47.net", "color-2048-seed-47.net",
        "color-3584-seed-47.net",  "color-512-seed-47.net",
        "color-6144-seed-47.net",  "color-7680-seed-47.net",
        "color-9216-seed-47.net",  "color-10240-seed-47.net",
        "color-11776-seed-47.net", "color-13312-seed-47.net",
        "color-14848-seed-47.net", "color-2560-seed-47.net",
        "color-4096-seed-47.net",  "color-5120-seed-47.net",
        "color-6656-seed-47.net",  "color-8192-seed-47.net",
        "color-9728-seed-47.net",  "color-10752-seed-47.net",
        "color-12288-seed-47.net", "color-13824-seed-47.net",
        "color-1536-seed-47.net",  "color-3072-seed-47.net",
        "color-4608-seed-47.net",  "color-5632-seed-47.net",
        "color-7168-seed-47.net",  "color-8704-seed-47.net",

        /*, "color-5120-seed-17.net",
        "color-5120-seed-19.net", "color-4608-seed-13.net",
        "color-4608-seed-17.net", "color-4608-seed-19.net",
        "color-4096-seed-13.net", "color-4096-seed-17.net",
        "color-4096-seed-19.net", "color-3584-seed-13.net",
        "color-3584-seed-17.net", "color-3584-seed-19.net",
        "color-3072-seed-13.net", "color-3072-seed-17.net",
        "color-3072-seed-19.net",*/
    };

    DOT::Solver solver;

    // Profiling algorithms
    std::array<double, 4> stat = {0, 0, 0, 0};
    std::array<double, 4> eapi = {0, 0, 0, 0};
    std::array<double, 4> eati = {0, 0, 0, 0};

    for (const auto &filename : Fs) {
      fprintf(stdout, "Step 1\n");
      fflush(stdout);

      solver.parseDimacs(base + SEP + filename);

      fprintf(stdout, "Step 2\n");
      fflush(stdout);

      // std::array<double, 4> times = solver.solveDimacs("bipartiteEapi");

      // eati[0] += times[0];
      // eati[1] += times[1];
      // eati[2] += times[2];
      // eati[3] += times[3];
      // stat[2] += times[0];

      //    stat[1] = solver.solveDimacs("colgen");

      stat[1] = solver.solveDimacs("bipartiteEapi");

      //    stat[1] = solver.solveDimacs("bipartiteEati");
      // eapi[0] += times[0];
      // eapi[1] += times[1];
      // eapi[2] += times[2];
      // eapi[3] += times[3];
      // stat[1] += times[0];

      //      stat[0] += solver.solveDimacs("netcplex");

      // stat[3] += solver.solverDimacs("colgen");
    }
  }

  if (false) {
    std::string SEP = "\\";

    std::string base =
        "C:\\Users\\Gualandi\\Google "
        "Drive\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";

    DOT::Histogram2D a;
    DOT::Histogram2D b;

    std::string f1 = "data32_1006.csv";
    std::string f2 = "data32_1007.csv";

    a.parse(base + "Shapes" + SEP + f1);
    b.parse(base + "Shapes" + SEP + f2);

    // a.parse("C:"
    //	"\\Users\\gualandi\\Documents\\GitHub\\dotlib\\nearby\\msvc\\Neraby"
    //	"\\x64\\Release\\test1.csv");
    // b.parse("C:"
    //	"\\Users\\gualandi\\Documents\\GitHub\\dotlib\\nearby\\msvc\\Neraby"
    //	"\\x64\\Release\\test2.csv");

    std::string msg = "Shapes";

    DOT::Solver solver;
    solver.init_coprimes(31, 32);
    // solver.tripartiteColgen(a, b);
    solver.bipartiteEapi(a, b, msg);
    solver.bipartiteEati(a, b, msg);
    solver.bipartiteCplex(a, b, msg);
  }

  if (false) {
    std::string SEP = "\\";

    std::string base =
        "C:\\Users\\Gualandi\\Google "
        "Drive\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";

    std::vector<std::string> dirs = {
        "ClassicImages",  //"Shapes"
        //"WhiteNoise", "CauchyDensity", "GRFmoderate","MicroscopyImages",
        //"GRFrough", "GRFsmooth", "LogGRF", "LogitGRF",,
    };

    std::vector<std::string> Fs = {
        "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
        "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv"};

    std::vector<std::string> Ss = {"128"};
    //, "64", "128", "256", "512"};
    std::array<double, 4> stat = {0, 0, 0, 0};
    std::array<double, 4> eapi = {0, 0, 0, 0};
    std::array<double, 4> eati = {0, 0, 0, 0};

    for (const auto &dtype : dirs) {
      for (const auto &S : Ss) {
        for (const auto &f11 : Fs) {
          for (const auto &f22 : Fs)
            if (f11 < f22) {
              DOT::Histogram2D a;
              DOT::Histogram2D b;

              std::string f1 = "data" + S + "_" + f11;
              std::string f2 = "data" + S + "_" + f22;

              a.parse(base + dtype + SEP + f1);
              b.parse(base + dtype + SEP + f2);

              std::string msg = dtype + " " + f11 + " " + f22;

              DOT::Solver solver;
              // solver.init_coprimes(31, atoi(S.c_str()));

              // for (int tau : {1,  5,  10, 15, 20, 25, 30, 35, 40, 45,
              //                50, 55, 60, 65, 70, 75, 80, 85, 95, 100})

              std::array<double, 4> times = solver.bipartiteEati(a, b, msg);
              eati[0] += times[0];
              eati[1] += times[1];
              eati[2] += times[2];
              eati[3] += times[3];
              stat[2] += times[0];

              times = solver.bipartiteEapi(a, b, msg);
              eapi[0] += times[0];
              eapi[1] += times[1];
              eapi[2] += times[2];
              eapi[3] += times[3];
              stat[1] += times[0];

              stat[0] += solver.bipartiteCplex(a, b, msg);

              stat[3] += solver.colgen(a, b, msg);
              fprintf(stdout, "\n");
            }
        }
      }
    }

    fprintf(stdout, "Cplex: %.4f\tEapi: %.4f\tEati: %.4f\tCG+eati: %.4f\n",
            stat[0], stat[1], stat[2], stat[3]);

    fprintf(stdout, "Eapi: %.4f\tpricing: %.4f\tbasis: %.4f\tduals: %.4f\n",
            eapi[0], eapi[1], eapi[2], eapi[3]);

    fprintf(stdout, "Eati: %.4f\tpricing: %.4f\tbasis: %.4f\tduals: %.4f\n",
            eati[0], eati[1], eati[2], eati[3]);
  }

  if (false) {
    size_t samples = n * n;
    DOT::Histogram2D a(n);
    DOT::Histogram2D b(n);

    vector<int64_t> tmp;
    int64_t cc = 0;
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++) {
        auto d = Uniform1N(gen);
        cc += d;
        tmp.push_back(d);
      }

    size_t idx = 0, m = tmp.size();
    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < n; j++) {
        a.add(i, j, tmp[idx]);
        b.add(i, j, tmp[m - idx - 1]);
        idx++;
      }

    PRINT("start solver %d %ld %ld\n", cc, a.balance(), b.balance());
    DOT::Solver solver;

    // double dist = solver.bipartite(a, b);

    // TODO: FIX double dist2 = solver.tripartite(a, b);
  }

  return EXIT_SUCCESS;
}
