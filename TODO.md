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


