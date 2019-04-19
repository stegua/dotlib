/* example1.c */


#include "emd.h"

#include <vector>
#include <chrono>

#include "DOT_BipartiteLemon.hpp"



namespace DOT {
// Ground distances
double dist(feature_t *F1, feature_t *F2) {
   int dX = F1->X - F2->X, dY = F1->Y - F2->Y;
   return static_cast<int64_t>(std::round(sqrt(dX*dX + dY * dY)* INTEGER_TOL));
   // NOTE: FOR COMPARISON MUST BE THE SAME AS_
   // #define DISTANCE_R2(x, y)  std::round(sqrt((x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]))*INTEGER_TOL)
}

// Use the Rubner code
void BipartiteEMD(const MeasureR2& Mu, const MeasureR2& Nu, int algo, const std::string& msg) {

   size_t m = Mu.size();
   size_t n = Nu.size();

   std::vector<double> w1(Mu.getWs());
   std::vector<double> w2(Nu.getWs());

   // Rubner signatures
   feature_t*   f1;
   feature_t*   f2;
   f1 = (feature_t*)malloc(m * sizeof(feature_t));
   f2 = (feature_t*)malloc(n * sizeof(feature_t));

   for (size_t i = 0; i < m; ++i) {
      auto p = Mu.getP(i);
      f1[i].X = p[0];
      f1[i].Y = p[1];
   }

   for (size_t j = 0; j < n; ++j) {
      auto p = Nu.getP(j);
      f2[j].X = p[0];
      f2[j].Y = p[1];
   }

   signature_t s1 = { m, f1, w1.data() },
               s2 = { n, f2, w2.data() };

   // Timinig output
   auto start = std::chrono::steady_clock::now();

   double       fobj;
   flow_t      *flow;
   int         flowSize;

   flow = (flow_t*)malloc(m * n * sizeof(flow_t));

   fobj = emd(&s1, &s2, dist, flow, &flowSize);

   free(flow);
   free(f1);
   free(f2);

   auto end = std::chrono::steady_clock::now();
   double elapsed = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

   fprintf(stdout, "%s %d %d Runtime %.6f Value %.f status %d\n",
           msg.c_str(), n, m, elapsed, fobj, 0);

   fflush(stdout);
}

}