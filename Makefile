solver:
	g++ -O3 -std=c++17 -march=native -DCINECA -o solver ./src/SolverCLI.cpp -I./include/ -I/marconi/home/userexternal/sgualand/solvers/lemon-1.3.1/ -lpthread -Wall
