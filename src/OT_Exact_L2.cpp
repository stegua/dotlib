/**
* @fileoverview Copyright (c) 2017-2018, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

#include "OT_Exact_L2.h""

/**
* @brief Solve standard Trasportation Problem on the complete bipartite graph
*        using Wassertein distance W_1^2 (order 1, with cost=L_2)
*/
real_t solve_exact_W_1_2(const histogram_t& h1, const histogram_t& h2) {
   auto logger = spd::get("console");

   using namespace lemon;

   real_t distance = std::numeric_limits<real_t>::max();

   size_t d = h1.size();
   size_t s = static_cast<size_t>(sqrt(d));

   auto ID = [&s](size_t x, size_t y) {
      return x * s + y;
   };

   // Time vars
   std::chrono::time_point<std::chrono::system_clock> start, end;
   // Start time.
   start = std::chrono::system_clock::now();

   // Build the graph for max flow
   Graph g;

   // add d nodes for each histrogam (d+1) source, (d+2) target
   std::vector<Graph::Node> nodes_A;
   nodes_A.reserve(d);
   // add first d source nodes for first partition
   for (size_t i = 0; i < d; ++i)
      nodes_A.emplace_back(g.addNode());

   std::vector<Graph::Node> nodes_B;
   nodes_B.reserve(d);
   // add first d source nodes for second partition
   for (size_t i = 0; i < d; ++i)
      nodes_B.emplace_back(g.addNode());

   // Add arcs for complete bipartite graph
   std::vector<Graph::Arc> arcs;
   arcs.reserve(d*d);
   std::vector<real_t> arcs_costs;
   arcs_costs.reserve(d*d);

   for (int i = 0; i < s; ++i)
      for (int j = 0; j < s; ++j) {
         for (int v = 0; v < s; ++v)
            for (int w = 0; w < s; ++w) {
               arcs.emplace_back(g.addArc(nodes_A[ID(i, j)], nodes_B[ID(v, w)]));
               arcs_costs.emplace_back(sqrt(pow(double(i - v), 2) + pow(double(j - w), 2)));
            }
      }

   // End time.
   end = std::chrono::system_clock::now();
   std::chrono::duration<double> inlineTimeElapsed = end - start;

   start = std::chrono::system_clock::now();

   NetworkSimplex<Graph, LimitValueType, real_t> cycle(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<real_t> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i) {
      b_i[nodes_A[i]] = +LimitValueType(h1[i]);
      b_i[nodes_B[i]] = -LimitValueType(h2[i]);
   }

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = cycle.INF;
      c_i[a] = std::trunc(MUL * arcs_costs[i]);
   }

   //set lower/upper bounds, cost
   cycle.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   NetworkSimplex<Graph, LimitValueType, real_t>::ProblemType ret = cycle.run();

   switch (ret) {
   case NetworkSimplex<Graph>::INFEASIBLE:
      logger->error("INFEASIBLE");
      break;
   case NetworkSimplex<Graph>::OPTIMAL:
      logger->info("OPTIMAL");
      break;
   case NetworkSimplex<Graph>::UNBOUNDED:
      logger->error("UNBOUNDED");
      break;
   }

   real_t sol_value = cycle.totalCost()/ MUL;

   end = std::chrono::system_clock::now();
   std::chrono::duration<double> run_time = end - start;

   logger->info("STE nodes {} arcs {} n {} N {} build_time {} cost {:.2f} run_time {}",
                countNodes(g), countArcs(g), h1.size(), s, inlineTimeElapsed.count(), sol_value, run_time.count());

   return sol_value;
}


/**
* @brief Solve standard Trasportation Problem on the complete bipartite graph
*        using Wassertein distance W_1^2 (order 1, with cost=L_2)
*/
real_t solve_exact_W_1_1(const histogram_t& h1, const histogram_t& h2) {
   auto logger = spd::get("console");

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
   std::vector<Graph::Node> nodes_A;
   nodes_A.reserve(d);
   // add first d source nodes for first partition
   for (size_t i = 0; i < d; ++i)
      nodes_A.emplace_back(g.addNode());

   std::vector<Graph::Node> nodes_B;
   nodes_B.reserve(d);
   // add first d source nodes for second partition
   for (size_t i = 0; i < d; ++i)
      nodes_B.emplace_back(g.addNode());

   // Add arcs for complete bipartite graph
   std::vector<Graph::Arc> arcs;
   arcs.reserve(d*d);
   std::vector<real_t> arcs_costs;
   arcs_costs.reserve(d*d);

   for (int i = 0; i < s; ++i)
      for (int j = 0; j < s; ++j) {
         for (int v = 0; v < s; ++v)
            for (int w = 0; w < s; ++w) {
               arcs.emplace_back(g.addArc(nodes_A[ID(i, j)], nodes_B[ID(v, w)]));
               arcs_costs.emplace_back(double(i - v) + double(j - w));
            }
      }

   logger->info("Input graph created with {} nodes and {} arcs", countNodes(g), countArcs(g));

   NetworkSimplex<Graph, LimitValueType, real_t> cycle(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<real_t> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i) {
      b_i[nodes_A[i]] = +LimitValueType(h1[i]);
      b_i[nodes_B[i]] = -LimitValueType(h2[i]);
   }

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = cycle.INF;
      c_i[a] = MUL*arcs_costs[i];
   }

   //set lower/upper bounds, cost
   cycle.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   NetworkSimplex<Graph, LimitValueType, real_t>::ProblemType ret = cycle.run();

   switch (ret) {
   case NetworkSimplex<Graph>::INFEASIBLE:
      logger->error("INFEASIBLE");
      break;
   case NetworkSimplex<Graph>::OPTIMAL:
      logger->info("OPTIMAL");
      break;
   case NetworkSimplex<Graph>::UNBOUNDED:
      logger->error("UNBOUNDED");
      break;
   }

   real_t sol_value = cycle.totalCost();

   return sol_value;
}


/**
* @brief Solve standard Trasportation Problem on the complete bipartite graph
*        using Wassertein distance W_1 ^ 2 (order 1, with cost = L_2)
*/
real_t solve_exact_W_1_8(const histogram_t& h1, const histogram_t& h2) {
   auto logger = spd::get("console");

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
   std::vector<Graph::Node> nodes_A;
   nodes_A.reserve(d);
   // add first d source nodes for first partition
   for (size_t i = 0; i < d; ++i)
      nodes_A.emplace_back(g.addNode());

   std::vector<Graph::Node> nodes_B;
   nodes_B.reserve(d);
   // add first d source nodes for second partition
   for (size_t i = 0; i < d; ++i)
      nodes_B.emplace_back(g.addNode());

   // Add arcs for complete bipartite graph
   std::vector<Graph::Arc> arcs;
   arcs.reserve(d*d);
   std::vector<real_t> arcs_costs;
   arcs_costs.reserve(d*d);

   for (int i = 0; i < s; ++i)
      for (int j = 0; j < s; ++j) {
         for (int v = 0; v < s; ++v)
            for (int w = 0; w < s; ++w) {
               arcs.emplace_back(g.addArc(nodes_A[ID(i, j)], nodes_B[ID(v, w)]));
               arcs_costs.emplace_back(std::max(double(i - v), double(j - w)));
            }
      }

   logger->info("Input graph created with {} nodes and {} arcs", countNodes(g), countArcs(g));

   NetworkSimplex<Graph, LimitValueType, real_t> cycle(g);

   // lower and upper bounds, cost
   ListDigraph::ArcMap<LimitValueType> l_i(g), u_i(g);
   ListDigraph::ArcMap<real_t> c_i(g);

   // FLow balance
   ListDigraph::NodeMap<LimitValueType> b_i(g);
   for (size_t i = 0; i < d; ++i) {
      b_i[nodes_A[i]] = +LimitValueType(h1[i]);
      b_i[nodes_B[i]] = -LimitValueType(h2[i]);
   }

   // Add all edges
   for (size_t i = 0, i_max = arcs.size(); i < i_max; ++i) {
      const auto& a = arcs[i];
      l_i[a] = 0;
      u_i[a] = cycle.INF;
      c_i[a] = std::trunc(10000 * 65536 * arcs_costs[i]) / (10000 * 65536);
   }

   //set lower/upper bounds, cost
   cycle.lowerMap(l_i).upperMap(u_i).costMap(c_i).supplyMap(b_i);

   NetworkSimplex<Graph, LimitValueType, real_t>::ProblemType ret = cycle.run();

   switch (ret) {
   case NetworkSimplex<Graph>::INFEASIBLE:
      logger->error("INFEASIBLE");
      break;
   case NetworkSimplex<Graph>::OPTIMAL:
      logger->info("OPTIMAL");
      break;
   case NetworkSimplex<Graph>::UNBOUNDED:
      logger->error("UNBOUNDED");
      break;
   }

   real_t sol_value = cycle.totalCost();

   return sol_value;
}