# Discrete Optimal Transport (DOT)

The main class to be developed are:

1. *Histogram2D*: 2D dimensional histogram type

2. *Config*: execution parameters and options


Functional perspective:

1. distance, plan <- compute_wdp(Hist2D a, Hist2D, config option)


2. distance, plan <- compute_wd1(Hist2D a, Hist2D, config option)

3. flow graph <- build_network_L1(Hist2D a, Hist2D)
4. flow graph <- build_network_Linf(Hist2D a, Hist2D)
5. flow graph <- build_network_L2(Hist2D a, Hist2D, int L) // if L=N-1, exact

6. distance, plan <- solve_min_cost_flow(flow graph)


7. Rispondere a questa domanda: https://stats.stackexchange.com/questions/404775/calculate-earth-movers-distance-for-two-grayscale-images

8. Check the input size: read by line, or all file at once

9. Compute Malhabois and Pearson distances to be used as ground distances in OT
cross	 check:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.mahalanobis.html

efficient implementation using cholesky:
https://stats.stackexchange.com/questions/65705/pairwise-mahalanobis-distances?rq=1

tutorial:
https://www.machinelearningplus.com/statistics/mahalanobis-distance/

linear equations:
https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/appendix-a-linear-solvers-basics/sparse-linear-systems/direct-method.html


