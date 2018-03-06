solver:
	g++ -O3 -std=c++17 -march=native -m64 -DCINECA -o solver ./src/SolverCLI.cpp -I./include/ -I./externs -I/marconi/home/userexternal/sgualand/solvers/lemon-1.3.1/ -lpthread -Wall

