/**
* @fileoverview Copyright (c) 2017-18, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once


#ifdef _WIN32
#include <amp.h>
#include <amp_math.h>

#include <ppl.h>
#endif

#include <chrono>
#include <random>

#include "DOT_BasicTypes.h"

extern "C" {
#include "gurobi_c.h"
#include "math.h"
}

#define POST(s, x) if (x) { fprintf(stdout,"%s: GRB error %d\n", s, x); exit(0); }

typedef std::pair<int, int>     Point2D;
typedef std::vector<Point2D>    Space2D;
typedef std::vector<double>     Marginal;

class Column {
 public:
   int     a;  // First point
   int     b;  // Second point
   double  c;  // Distance

   Column() : a(-1), b(-1), c(-1) {}
   Column(int _a, int _b, double _c) : a(_a), b(_b), c(_c) {}
};

typedef std::vector<Column>                Columns;


namespace DOT {

Space2D uniformSpace2D(size_t n, int seed) {
   std::random_device rd;  //Will be used to obtain a seed for the random number engine
   std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
   std::uniform_int_distribution<> Uniform0N(0, n);

   Space2D r;
   r.reserve(n*n);
   for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
         r.emplace_back(i,j);

   return r;
}

Space2D randomSpace2D(size_t n, int seed) {
   std::random_device rd;  //Will be used to obtain a seed for the random number engine
   std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
   std::uniform_int_distribution<> Uniform0N(0, n);

   Space2D r;
   r.reserve(n);
   for (size_t i = 0; i < n; ++i)
      r.emplace_back(Uniform0N(gen), Uniform0N(gen));
   return r;
}

Marginal randomMarginal(size_t n, int max_value, int seed) {
   std::random_device rd;  //Will be used to obtain a seed for the random number engine
   std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
   std::uniform_int_distribution<> Uniform1N(1, max_value);

   Marginal pi;
   pi.reserve(n);
   double tot = 0;
   for (size_t i = 0; i < n; ++i) {
      double val = Uniform1N(gen);
      tot += val;
      pi.emplace_back(val);
   }
   for (auto& v : pi)
      v = v / tot;

   return pi;
}


double distance(const Point2D& a, const Point2D& b) {
   return sqrt(pow((a.first - b.first), 2) + pow((a.second - b.second), 2));
}

Columns HeuristicColumns(const Space2D& X, const Space2D& Y, double threshold) {
   Columns E;

   // Assumption X e Y are different spaces
   E.reserve(X.size()*Y.size());

   for (int i = 0, i_max = X.size(); i < i_max; ++i)
      for (int j = 0, j_max = Y.size(); j < j_max; ++j) {
         double c = distance(X[i], Y[j]);
         if (c <= threshold)
            E.emplace_back(i, j, c);
      }

   return E;
}

std::pair<double, Column> solvePricing(const Space2D& X,
                                       const std::vector<double>& DualsPI,
                                       const Space2D& Y,
                                       const std::vector<double>& DualsMU) {
   int a = -1;
   int b = -1;
   double c = -1;
   double viol = -1;

   for (int i = 0, i_max = X.size(); i < i_max; ++i) {
      for (int j = 0, j_max = Y.size(); j < j_max; ++j) {
         double d = distance(X[i], Y[j]);
         double v1 = DualsPI[i] + DualsMU[j] - d;
         if (v1 > viol) {
            a = i;
            b = j;
            c = d;
            viol = v1;
         }
      }
   }

   return std::make_pair(viol, Column(a,b,c));
}

std::pair<double, Columns> solvePricing2(
   const Space2D& X,
   const std::vector<double>& DualsPI,
   const Space2D& Y,
   const std::vector<double>& DualsMU) {
   Columns cols;

   double max_viol = -1;
   for (int i = 0, i_max = X.size(); i < i_max; ++i) {
      int a = -1;
      int b = -1;
      double c = -1;
      double viol = -1;
      for (int j = 0, j_max = Y.size(); j < j_max; ++j) {
         double d = distance(X[i], Y[j]);
         double v1 = DualsPI[i] + DualsMU[j] - d;
         if (v1 > viol) {
            a = i;
            b = j;
            c = d;
            viol = v1;
            if (v1 > max_viol)
               max_viol = v1;
         }
      }
      if (viol > 0.001)
         cols.emplace_back(a, b, c);
   }

   return std::make_pair(max_viol, cols);
}

void solvePricing3(
   const Space2D& X,
   const std::vector<double>& DualsPI,
   const Space2D& Y,
   const std::vector<double>& DualsMU,
   Columns& out) {
   concurrency::parallel_for(0,(int) X.size(),
   [&](int i) {
      int b = -1;
      double c = -1;
      double viol = -1;
      // double or float?
      for (int j = 0, j_max = Y.size(); j < j_max; ++j) {
         double v1 = DualsPI[i] + DualsMU[j];
         if (v1 > 0.99) {
            double d = distance(X[i], Y[j]);
            v1 = v1 - d;
            if (v1 > viol) {
               b = j;
               c = d;
               viol = v1;
            }
         }
         if (viol > 0.001) {
            out[i].b = b;
            out[i].c = c;
         }
      }
   });
}

void solvePricing4(
   int N,
   const concurrency::array<Point2D>& xv,
   const std::vector<double>& DualsPI,
   const concurrency::array<Point2D>& yv,
   const std::vector<double>& DualsMU,
   Columns& cols) {

   concurrency::array<double> pi(DualsPI.size(), DualsPI.begin(), DualsPI.end());
   concurrency::array<double> mu(DualsMU.size(), DualsMU.begin(), DualsMU.end());

   concurrency::array_view<Column> cv(cols.size(), cols);

   concurrency::parallel_for_each(cv.extent, [=, &xv, &yv, &pi, &mu](concurrency::index<1> idx) restrict(amp) {
      int b = -1;
      float c = -1;
      float viol = -1;
      for (int j = 0; j < N; ++j) {
         float v1 = pi[idx] + mu[j];
         if (v1 > 0.99) {
            float d = concurrency::fast_math::sqrt(concurrency::fast_math::pow(xv[idx].first - yv[j].first, 2) + concurrency::fast_math::pow(xv[idx].second - yv[j].second, 2));
            v1 = v1 - d;
            if (v1 > viol) {
               b = j;
               c = d;
               viol = v1;
            }
         }
      }
      if (viol > 0.001) {
         cv[idx].b = b;
         cv[idx].c = c;
      }
   });
}

void ColumnGeneration(const Space2D& X, const Marginal& PI,
                      const Space2D& Y, const Marginal& MU) {
   auto start = std::chrono::steady_clock::now();
   double elapsed;

   // Condition on input data
   assert(X.size() == PI.size() && Y.size() == MU.size());

   size_t n = PI.size();
   size_t m = MU.size();

   // Pointers for Gurobi
   GRBenv   *env_master;
   GRBmodel *master;
   // Arrays for passing data to Gurobi
   int     ind[2];
   double  val[2] = { 1, 1 };

   // Build master problem
   GRBloadenv(&env_master, NULL);
   GRBsetintparam(env_master, GRB_INT_PAR_OUTPUTFLAG, 1);
   GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
   GRBsetintparam(env_master, GRB_INT_PAR_PRESOLVE, 0);
   //GRBsetintparam(env_master, GRB_INT_PAR_PRECRUSH, 1);  // TO BE SET TO 1 SINCE WE USE USER CUTS
   GRBsetintparam(env_master, GRB_INT_PAR_THREADS, 1);
   GRBsetdblparam(env_master, GRB_DBL_PAR_TIMELIMIT, 60); // 2 minuti

   GRBnewmodel(env_master, &master, "KantorovichCG", 0, NULL, NULL, NULL, NULL, NULL);
   GRBsetintattr(master, "ModelSense", GRB_MINIMIZE);

   // Empty constraints
   for (const auto& pi: PI)
      POST("add supply constraint", GRBaddconstr(master, 0, ind, val, GRB_GREATER_EQUAL, pi, ""));

   for (const auto& mu : MU)
      POST("add demand constraint", GRBaddconstr(master, 0, ind, val, GRB_LESS_EQUAL, mu, ""));

   GRBupdatemodel(master);

   // Heuristic set of columns
   int n_cols = 0;
   Columns E = HeuristicColumns(X, Y, 2);
   for (const Column& col : E) {
      n_cols++;
      ind[0] = col.a;
      ind[1] = n + col.b;
      POST("add column constraint", GRBaddvar(master, 2, ind, val, col.c, 0, GRB_INFINITY, GRB_CONTINUOUS, ""));
   }

   // Timinig output
   auto end = std::chrono::steady_clock::now();
   elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
   printf("Basic model up: %.4f\n", elapsed);

   // Solve the problem
   start = std::chrono::steady_clock::now();

   double   obj;
   int      status;
   int it = 0;
   std::vector<double> DualsPI(n, 0);
   std::vector<double> DualsMU(n, 0);

   concurrency::array<Point2D> xv(X.size(), X.begin(), X.end());
   concurrency::array<Point2D> yv(Y.size(), Y.begin(), Y.end());

   Columns cols(n);
   for (int i = 0; i < n; ++i) {
      cols[i].a = i;
      cols[i].c = -1;
   }

   while (it++ < 100000) {
      GRBoptimize(master);

      GRBgetintattr(master, "Status", &status);
      if (status == GRB_UNBOUNDED || status == GRB_INFEASIBLE) {
         fprintf(stdout, "Unbounded or Infeasible\n");
         goto QUIT;
      }

      // Take the current LP decision vector
      POST("get obj", GRBgetdblattr(master, "ObjVal", &obj));

      end = std::chrono::steady_clock::now();
      elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
      printf("it %d: Time %.3f - Value:%.6f - NumVars: %d\n", it, elapsed, obj, n_cols);

      // Take the dual values
      GRBgetdblattrarray(master, "Pi", 0, n, &DualsPI[0]);
      GRBgetdblattrarray(master, "Pi", n, m, &DualsMU[0]);

      //solvePricing3(X, DualsPI, Y, DualsMU, cols);
      solvePricing4(Y.size(), xv, DualsPI, yv, DualsMU, cols);

      // Add new column
      bool stop = true;
      for (Column& col : cols)
         if (col.c > 0) {
            stop = false;
            n_cols++;
            ind[0] = col.a;
            ind[1] = n + col.b;
            POST("add column constraint", GRBaddvar(master, 2, ind, val, col.c, 0, GRB_INFINITY, GRB_CONTINUOUS, ""));
            col.c = -1; // reset cost value
         }
      if (stop) break;
   }
   end = std::chrono::steady_clock::now();
   elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
   printf("Final Time %.3f - Value: %.3f - NumVars: %d\n", elapsed, obj, n_cols);

QUIT:
   POST("free", GRBfreemodel(master));
   GRBfreeenv(env_master);
}

} /// END Namespace DOT