/**
* @fileoverview Copyright (c) 2017, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#pragma once

#include "OT_BasicTypes.h"

/**
* @brief Compute Earth Moving Distance (EMD) between pair of images
*        spacing razs at "constant" angle distance, with given degree ration
*        Example: with degree=45, it puts 8 rays for vertex (at most)
*/
real_t apx_L2_1(const histogram_t& h1, const histogram_t& h2, double degree = 1) {
   using namespace lemon;

   size_t d = h1.size();
   size_t s = static_cast<size_t>(sqrt(d));

   auto ID = [&s](size_t x, size_t y) {
      return x * s + y;
   };

   auto cart2pol = [&s](double x, double y) {
      return atan2(y, x);
   };

   // Build the graph for max flow
   Graph g;

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<Graph::Node> nodes;
   std::vector<std::pair<int, int>> points;
   points.reserve(d);
   for (size_t i = 0; i < s; ++i)
      for (size_t j = 0; j < s; ++j) {
         points.emplace_back(std::make_pair(i, j));
         nodes.emplace_back(g.addNode());
      }
   assert(points.size() == d);

   std::vector<Graph::Arc> arcs;
   std::vector<real_t> arcs_costs;

   pair_set A;

   // First iterate for every node of the graph
   size_t max_w = s;
   for (size_t l = 0; l < d; ++l) {
      int_set D;
      int v = points[l].first;
      int w = points[l].second;
      // fprintf(stdout, "point %d %d\n", v, w);
      for (int i = v; i < v + max_w; ++i)
         for (int j = w; j < w + max_w; ++j) {
            if (!(i == v && j == w)) {
               if (i < s && j < s) {
                  int f = int(lround(double(180.0 / degree / M_PI * cart2pol(i - v, j - w))) * degree);
                  if (D.find(f) == D.end()) {
                     D.insert(f);

                     arcs.emplace_back(g.addArc(nodes[ID(v, w)], nodes[ID(i, j)]));
                     A.insert(std::make_pair(ID(v, w), ID(i, j)));
                     arcs_costs.emplace_back(sqrt(pow(i - v, 2) + pow(j - w, 2)));
//                     fprintf(stdout, "(%d,%d) -> (%d,%d)#%.3f\t", v, w, i, j, arcs_costs.back());

                     arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(v, w)]));
                     A.insert(std::make_pair(ID(i, j), ID(v, w)));
                     arcs_costs.emplace_back(sqrt(pow(i - v, 2) + pow(j - w, 2)));
                     //fprintf(stdout, "(%d,%d) -> (%d,%d)#%.3f\t", i, j, v, w, arcs_costs.back());
                  }
               }
            }
         }

      // up down
      for (int i = 1; i < max_w; ++i)
         for (int j = w + 1; j < w + max_w; ++j) {
            int aa = v;
            int bb = w;
            int cc = v - i;
            int dd = j;
            int f = int(lround(double(180.0 / degree / M_PI * cart2pol(aa - cc, dd - bb))) * degree - 90);
            //fprintf(stdout, "%d %d %d %d angle: %d\n", aa, bb, cc, dd, f);
            if (!(aa == cc && bb == dd)) {
               if (cc >= 0 && dd < s) {
                  if (D.find(f) == D.end()) {
                     D.insert(f);
                     arcs.emplace_back(g.addArc(nodes[ID(aa, bb)], nodes[ID(cc, dd)]));
                     A.insert(std::make_pair(ID(aa, bb), ID(cc, dd)));
                     arcs_costs.emplace_back(sqrt(pow(aa - cc, 2) + pow(dd - bb, 2)));
                     //       fprintf(stdout, "(%d,%d) -> (%d,%d)#%.3f\t", aa, bb, cc, dd, arcs_costs.back());

                     arcs.emplace_back(g.addArc(nodes[ID(cc, dd)], nodes[ID(aa, bb)]));
                     A.insert(std::make_pair(ID(cc, dd), ID(aa, bb)));
                     arcs_costs.emplace_back(sqrt(pow(aa - cc, 2) + pow(dd - bb, 2)));
                     //     fprintf(stdout, "(%d,%d) -> (%d,%d)#%.3f\t", cc, dd, aa, bb, arcs_costs.back());
                  }
               }
            }
         }

      //if (l > s)
      //   exit(0);
      //fprintf(stdout, "\n");
   }


   fprintf(stdout, "Input graph created with %d nodes and %d arcs - check: %d\n", countNodes(g), countArcs(g), A.size());
   for (auto d : arcs_costs)
      if (d < 1)
         fprintf(stdout, "%.3f ", d);
//exit(0);
   NetworkSimplex<Graph, LimitValueType, real_t> cycle(g);

// lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<real_t> c_i(g);

// FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i)
      b_i[nodes[i]] = int(h1[i] - h2[i]);

// Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = cycle.INF;
      c_i[a] = std::trunc(65536 * arcs_costs[i]) / 65536;
   }

//set lower/upper bounds, cost
   cycle.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   NetworkSimplex<Graph, LimitValueType, real_t>::ProblemType ret = cycle.run();

   switch (ret) {
   case NetworkSimplex<Graph>::INFEASIBLE:
      std::cerr << "INFEASIBLE" << std::endl;
      break;
   case NetworkSimplex<Graph>::OPTIMAL:
      std::cerr << "OPTIMAL" << std::endl;
      break;
   case NetworkSimplex<Graph>::UNBOUNDED:
      std::cerr << "UNBOUNDED" << std::endl;
      break;
   }

   real_t sol_value = cycle.totalCost();

   return sol_value;


}

/**
* @brief Compute Earth Moving Distance (EMD) between pair of images
*/
real_t apx_L2_2(const histogram_t& h1, const histogram_t& h2) {
   using namespace lemon;

   real_t distance = std::numeric_limits<real_t>::max();

   size_t d = h1.size();
   size_t s = static_cast<size_t>(sqrt(d));

   auto ID = [&s](size_t x, size_t y) {
      return x * s + y;
   };

   // Build the graph for max flow
   Graph g;

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<Graph::Node> nodes;
   // add first d source nodes
   for (size_t i = 0; i < d; ++i)
      nodes.emplace_back(g.addNode());

   std::vector<Graph::Arc> arcs;
   std::vector<real_t> arcs_costs;
   arcs.reserve(4 * d);
   for (size_t i = 0; i < s; ++i)
      for (size_t j = 0; j < s - 1; ++j) {
         arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i, j + 1)]));
         arcs.emplace_back(g.addArc(nodes[ID(i, j + 1)], nodes[ID(i, j)]));

         arcs_costs.emplace_back(1.0);
         arcs_costs.emplace_back(1.0);
      }

   for (size_t i = 0; i < s - 1; ++i)
      for (size_t j = 0; j < s; ++j) {
         arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 1, j)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 1, j)], nodes[ID(i, j)]));

         arcs_costs.emplace_back(1.0);
         arcs_costs.emplace_back(1.0);
      }

   for (size_t i = 0; i < s - 1; ++i)
      for (size_t j = 0; j < s - 1; ++j) {
         arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 1, j + 1)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 1, j + 1)], nodes[ID(i, j)]));
         arcs.emplace_back(g.addArc(nodes[ID(i, j + 1)], nodes[ID(i + 1, j)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 1, j)], nodes[ID(i, j + 1)]));

         arcs_costs.emplace_back(sqrt(2.0));
         arcs_costs.emplace_back(sqrt(2.0));
         arcs_costs.emplace_back(sqrt(2.0));
         arcs_costs.emplace_back(sqrt(2.0));
      }

   for (size_t i = 0; i < s - 2; ++i)
      for (size_t j = 0; j < s - 2; ++j) {
         arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 2, j + 1)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 2, j + 1)], nodes[ID(i, j)]));

         arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 1, j + 2)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 1, j + 2)], nodes[ID(i, j)]));

         arcs.emplace_back(g.addArc(nodes[ID(i, j + 2)], nodes[ID(i + 1, j)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 1, j)], nodes[ID(i, j + 2)]));

         arcs.emplace_back(g.addArc(nodes[ID(i, j + 2)], nodes[ID(i + 2, j + 1)]));
         arcs.emplace_back(g.addArc(nodes[ID(i + 2, j + 1)], nodes[ID(i, j + 2)]));

         arcs_costs.emplace_back(sqrt(5.0));
         arcs_costs.emplace_back(sqrt(5.0));
         arcs_costs.emplace_back(sqrt(5.0));
         arcs_costs.emplace_back(sqrt(5.0));
         arcs_costs.emplace_back(sqrt(5.0));
         arcs_costs.emplace_back(sqrt(5.0));
         arcs_costs.emplace_back(sqrt(5.0));
         arcs_costs.emplace_back(sqrt(5.0));
      }
   //for (size_t i = 0; i < s - 3; ++i)
   //   for (size_t j = 0; j < s - 3; ++j) {
   //      arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 3, j + 1)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 3, j + 1)], nodes[ID(i, j)]));

   //      arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 1, j + 3)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 1, j + 3)], nodes[ID(i, j)]));

   //      arcs.emplace_back(g.addArc(nodes[ID(i, j + 3)], nodes[ID(i + 1, j)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 1, j)], nodes[ID(i, j + 3)]));

   //      arcs.emplace_back(g.addArc(nodes[ID(i, j + 3)], nodes[ID(i + 3, j + 2)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 3, j + 2)], nodes[ID(i, j + 3)]));

   //      arcs_costs.emplace_back(sqrt(10.0));
   //      arcs_costs.emplace_back(sqrt(10.0));
   //      arcs_costs.emplace_back(sqrt(10.0));
   //      arcs_costs.emplace_back(sqrt(10.0));
   //      arcs_costs.emplace_back(sqrt(10.0));
   //      arcs_costs.emplace_back(sqrt(10.0));
   //      arcs_costs.emplace_back(sqrt(10.0));
   //      arcs_costs.emplace_back(sqrt(10.0));
   //   }

   //for (size_t i = 0; i < s - 4; ++i)
   //   for (size_t j = 0; j < s - 4; ++j) {
   //      arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 4, j + 1)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 4, j + 1)], nodes[ID(i, j)]));

   //      arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 1, j + 4)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 1, j + 4)], nodes[ID(i, j)]));

   //      arcs.emplace_back(g.addArc(nodes[ID(i, j + 4)], nodes[ID(i + 1, j)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 1, j)], nodes[ID(i, j + 4)]));

   //      arcs.emplace_back(g.addArc(nodes[ID(i, j + 4)], nodes[ID(i + 4, j + 3)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 4, j + 3)], nodes[ID(i, j + 4)]));

   //      arcs_costs.emplace_back(sqrt(17.0));
   //      arcs_costs.emplace_back(sqrt(17.0));
   //      arcs_costs.emplace_back(sqrt(17.0));
   //      arcs_costs.emplace_back(sqrt(17.0));
   //      arcs_costs.emplace_back(sqrt(17.0));
   //      arcs_costs.emplace_back(sqrt(17.0));
   //      arcs_costs.emplace_back(sqrt(17.0));
   //      arcs_costs.emplace_back(sqrt(17.0));
   //   }

   //for (size_t i = 0; i < s - 4; ++i)
   //   for (size_t j = 0; j < s - 4; ++j) {
   //      arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 4, j + 3)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 4, j + 3)], nodes[ID(i, j)]));

   //      arcs.emplace_back(g.addArc(nodes[ID(i, j)], nodes[ID(i + 3, j + 4)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 3, j + 4)], nodes[ID(i, j)]));

   //      arcs.emplace_back(g.addArc(nodes[ID(i, j + 4)], nodes[ID(i + 3, j)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 3, j)], nodes[ID(i, j + 4)]));

   //      arcs.emplace_back(g.addArc(nodes[ID(i, j + 4)], nodes[ID(i + 4, j + 1)]));
   //      arcs.emplace_back(g.addArc(nodes[ID(i + 4, j + 1)], nodes[ID(i, j + 4)]));

   //      arcs_costs.emplace_back(sqrt(25.0));
   //      arcs_costs.emplace_back(sqrt(25.0));
   //      arcs_costs.emplace_back(sqrt(25.0));
   //      arcs_costs.emplace_back(sqrt(25.0));
   //      arcs_costs.emplace_back(sqrt(25.0));
   //      arcs_costs.emplace_back(sqrt(25.0));
   //      arcs_costs.emplace_back(sqrt(25.0));
   //      arcs_costs.emplace_back(sqrt(25.0));
   //   }



   fprintf(stdout, "Input graph created with %d nodes and %d arcs\n",
           countNodes(g), countArcs(g));

   NetworkSimplex<Graph, LimitValueType, real_t> cycle(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<real_t> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i)
      b_i[nodes[i]] = int(h1[i] - h2[i]);

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = cycle.INF;
      c_i[a] = std::trunc(65536 * arcs_costs[i]) / 65536;
   }

   //set lower/upper bounds, cost
   cycle.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   NetworkSimplex<Graph, LimitValueType, real_t>::ProblemType ret = cycle.run();

   switch (ret) {
   case NetworkSimplex<Graph>::INFEASIBLE:
      std::cerr << "INFEASIBLE" << std::endl;
      break;
   case NetworkSimplex<Graph>::OPTIMAL:
      std::cerr << "OPTIMAL" << std::endl;
      break;
   case NetworkSimplex<Graph>::UNBOUNDED:
      std::cerr << "UNBOUNDED" << std::endl;
      break;
   }

   real_t sol_value = cycle.totalCost();

   return sol_value;
}