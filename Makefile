#----------------------------------------------------------
# DOTLib v0.4.0: Discrete Optimal Transport library
#
# @fileoverview Copyright (c) 2017-20, Stefano Gualandi,
#               via Ferrata, 1, I-27100, Pavia, Italy
#
# @author stefano.gualandi@gmail.com (Stefano Gualandi)
#----------------------------------------------------------

OPTFLAG = -O2 -ffast-math -DNDEBUG -Wall -std=c++17 -DLINUX
LDFLAGS = -O2 -DNDEBUG -pthread -std=c++17 -lm

COMPILER    = g++ ${OPTFLAG}
LINKER      = g++ ${LDFLAGS}

solver::
	${COMPILER} -o solver.o -c ./src/SolverW1L2.cpp -I./include/ -I./externs
	${LINKER} -o solver ./solver.o

convex::
	${COMPILER} -DCONVEXHULL -o solver.o -c ./src/SolverW1L2.cpp -I./include/ -I./externs
	${LINKER} -o solver ./solver.o

test::
	./solver -f data/test1.csv -g data/test3.csv -L 1 
	./solver -f data/test1.csv -g data/test3.csv -L 2
	./solver -f data/test1.csv -g data/test3.csv -L 3

dotmark32::
	./solver -f data/data32_1001.csv -g data/data32_1002.csv -L 2
	./solver -f data/data32_1001.csv -g data/data32_1002.csv -L 3
	./solver -f data/data32_1001.csv -g data/data32_1002.csv -L 5
	./solver -f data/data32_1001.csv -g data/data32_1002.csv -L 10
	./solver -f data/data32_1001.csv -g data/data32_1002.csv -L 31

dotmark128::
	./solver -f data/data128_1001.csv -g data/data128_1002.csv -L 2
	./solver -f data/data128_1001.csv -g data/data128_1002.csv -L 3
	./solver -f data/data128_1001.csv -g data/data128_1002.csv -L 5
	./solver -f data/data128_1001.csv -g data/data128_1002.csv -L 10
	
clean::
	rm ./solver ./solver.o
