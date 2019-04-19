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

// Gurobi LP library
extern "C" {
#include "gurobi_c.h"
#include "math.h"
}

#include "DOT_MeasureR2.hpp"

#define POST(s, x) if (x) { fprintf(stdout,"%s: GRB error %d\n", s, x); exit(0); }


namespace DOT {

// Compute Kantorovich-Wasserstein distance between two measures defined in R2 (bipartite version
void BipartiteLP(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg) {
   auto start = std::chrono::steady_clock::now();
   double elapsed;

   int m = (int)Mu.size();
   int n = (int)Nu.size();

   // Pointers for Gurobi
   GRBenv   *env_master;
   GRBmodel *master;

   // Build master problem
   GRBloadenv(&env_master, NULL);

   if (algo == 0)
      GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_PRIMAL);
   if (algo == 1)
      GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
   if (algo == 2) {
      GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER);
      GRBsetintparam(env_master, GRB_INT_PAR_CROSSOVER, 0);
   }

   GRBsetintparam(env_master, GRB_INT_PAR_OUTPUTFLAG, 0);
   GRBsetdblparam(env_master, GRB_DBL_PAR_FEASIBILITYTOL, 1e-09);

   // Arrays for passing data to Gurobi
   vector<double>  obj(m*n, 0);
   for (size_t i = 0; i < m; ++i)
      for (size_t j = 0; j < n; ++j)
         obj[i*m+j] = (double)DISTANCE_R2(Mu.getP(i), Nu.getP(j));


   GRBnewmodel(env_master, &master, "KWD", m*n, &obj[0], NULL, NULL, NULL, NULL);
   GRBsetintattr(master, "ModelSense", GRB_MINIMIZE);

   GRBupdatemodel(master);

   // Full set of constraints
   vector<int>     ind(m, 0);
   vector<double>  val(m, 1);

   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
         ind[j] = m*i + j;
      }
      POST("add cut constraint", GRBaddconstr(master, n, &ind[0], &val[0], GRB_EQUAL, (double)Mu.getW(i), NULL));
   }

   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
         ind[j] = m * j + i;
      }
      POST("add cut constraint", GRBaddconstr(master, n, &ind[0], &val[0], GRB_EQUAL, (double)Nu.getW(i), NULL));
   }

   GRBupdatemodel(master);

   // Timinig output
   auto end = std::chrono::steady_clock::now();
   elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   // Solve the problem
   start = std::chrono::steady_clock::now();

   double  fobj = 0.0;
   double  gub_tot = 0.0;
   int     status;

   GRBoptimize(master);

   GRBgetintattr(master, "Status", &status);

   if (status == GRB_UNBOUNDED || status == GRB_INFEASIBLE) {
      fprintf(stdout, "Unbounded or Infeasible\n");
      goto QUIT;
   }

   // Take the current LP decision vector
   POST("get obj", GRBgetdblattr(master, "ObjVal", &fobj));
   POST("get runtime", GRBgetdblattr(master, "Runtime", &gub_tot));

   end = std::chrono::steady_clock::now();
   elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   fprintf(stdout, "%s %d %d Runtime %.6f Value %.f status %d\n",
           msg.c_str(), n, m, elapsed, fobj, status);

   fflush(stdout);

QUIT:
   POST("free", GRBfreemodel(master));
   GRBfreeenv(env_master);
}

// Compute Kantorovich-Wasserstein distance between two measures defined in R2 (bipartite version
void BipartiteLP_Dual(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg) {
   auto start = std::chrono::steady_clock::now();
   double elapsed;

   int m = (int)Mu.size();
   int n = (int)Nu.size();

   // Pointers for Gurobi
   GRBenv   *env_master;
   GRBmodel *master;

   // Build master problem
   GRBloadenv(&env_master, NULL);

   if (algo == 0)
      GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_PRIMAL);
   if (algo == 1)
      GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
   if (algo == 2) {
      GRBsetintparam(env_master, GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER);
      GRBsetintparam(env_master, GRB_INT_PAR_CROSSOVER, 0);
   }

   GRBsetintparam(env_master, GRB_INT_PAR_OUTPUTFLAG, 0);
   GRBsetdblparam(env_master, GRB_DBL_PAR_FEASIBILITYTOL, 1e-09);

   // Arrays for passing data to Gurobi
   vector<double>  obj(n + m, 0);
   vector<double>  vlb(n + m, 0);
   vector<double>  vub(n + m, 0);
   for (size_t i = 0; i < m; ++i) {
      obj[i] = Mu.getW(i);
      vlb[i] = -GRB_INFINITY;
      vub[i] = +GRB_INFINITY;
   }

   for (size_t i = 0; i < n; ++i) {
      obj[m + i] = Nu.getW(i);
      vlb[m + i] = -GRB_INFINITY;
      vub[m + i] = +GRB_INFINITY;
   }

   GRBnewmodel(env_master, &master, "KWD", int(m + n), &obj[0], &vlb[0], &vub[0], NULL, NULL);
   GRBsetintattr(master, "ModelSense", GRB_MAXIMIZE);

   GRBupdatemodel(master);

   // Full set of constraints
   int     ind[2] = { 0, 0 };
   double  val[2] = { 1, 1 };
   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
         ind[0] = i;
         ind[1] = m + j;
         POST("add cut constraint", GRBaddconstr(master, 2, ind, val, GRB_LESS_EQUAL, (double)DISTANCE_R2(Mu.getP(i), Nu.getP(j)), NULL));
      }
   }

   GRBupdatemodel(master);

   // Timinig output
   auto end = std::chrono::steady_clock::now();
   elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   // Solve the problem
   start = std::chrono::steady_clock::now();

   double  fobj = 0.0;
   double  gub_tot = 0.0;
   int     status;

   GRBoptimize(master);

   GRBgetintattr(master, "Status", &status);

   if (status == GRB_UNBOUNDED || status == GRB_INFEASIBLE) {
      fprintf(stdout, "Unbounded or Infeasible\n");
      goto QUIT;
   }

   // Take the current LP decision vector
   POST("get obj", GRBgetdblattr(master, "ObjVal", &fobj));
   POST("get runtime", GRBgetdblattr(master, "Runtime", &gub_tot));

   end = std::chrono::steady_clock::now();
   elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   fprintf(stdout, "%s %d %d Runtime %.6f Value %.f status %d\n",
           msg.c_str(), n, m, elapsed, fobj, status);

   fflush(stdout);

QUIT:
   POST("free", GRBfreemodel(master));
   GRBfreeenv(env_master);
}

} // End namespace