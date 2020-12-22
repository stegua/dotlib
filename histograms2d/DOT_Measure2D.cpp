/**
 * @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

 #include "DOT_Measure1D.h"

namespace DOT {
// // Solve Optimal Transport between two 1D measures
// std::pair<double, Plan> Wasserstein2_1D(const Measure1D& mu, const Measure1D& nu) {
//     Plan pi;        // Trasport Plan
//     double z = 0;   // Optimal Value

//     size_t m = mu.size();
//     size_t n = nu.size();

//     double a = mu.getWeight(0);    
//     double b = nu.getWeight(0);

//     for (size_t i = 0, j = 0; i < m && j < n; ) {
//       if (a == b) {
//             z += a * pow(mu.getPoint(i) - nu.getPoint(j), 2);
//             pi.add(i, j, a);
//             a = mu.getWeight(++i);
//             b = nu.getWeight(++j);
//       } else {
//             double f = a < b ? a : b;
//             z += f * pow(mu.getPoint(i) - nu.getPoint(j), 2);
//             pi.add(i, j, f);
//             if (a > b) {
//                 a = a - b;
//                 b = nu.getWeight(++j);
//             }
//             else {
//                 b = b - a;
//                 a = mu.getWeight(++i);
//         }
//     }
//   }

//     return std::make_pair(z, pi);
//   }

}