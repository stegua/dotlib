/*
 * @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
 *               via Ferrata, 5, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#include <random>

#include "DOT_Histogram2D.h"

int main(int argc, char* argv[]) {
	int n = 8;

	if (argc > 1)
		n = atoi(argv[1]);

	int seed = 13;

	std::random_device
		rd; // Will be used to obtain a seed for the random number engine
	std::mt19937 gen(seed); // Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> Uniform01(0, 1);
	std::uniform_int_distribution<> Uniform1N(1, n);

	if (true) {
		DOT::Histogram2D a;
		DOT::Histogram2D b;

		a.parse("C:\\Users\\Gualandi\\Google "
			"Drive\\Ricerca\\DOTA\\data\\DOTmark_1."
			"0\\Data\\ClassicImages\\data128_1001.csv");
		b.parse("C:\\Users\\Gualandi\\Google "
			"Drive\\Ricerca\\DOTA\\data\\DOTmark_1."
			"0\\Data\\ClassicImages\\data128_1002.csv");

		//std::string pp = "C:\\Users\\Gualandi\\Documents\\GitHub\\dotlib\\nearby\\msvc\\Neraby\\x64\\Debug\\";
		//a.parse(pp + "test1.csv");
		//b.parse(pp + "test2.csv");


		PRINT("start solver %ld %ld\n", a.balance(), b.balance());
		DOT::Solver solver;

		//double dist0 = solver.bipartite(a, b);

		// TODO: FIX double dist1 = solver.tripartite(a, b);

		//solver.phaseOne(a, b);

		//solver.phaseTwo(a, b);

		solver.colgen(a, b);
	}

	if (false) {
		size_t samples = n * n;
		DOT::Histogram2D a(n);
		DOT::Histogram2D b(n);

		vector<int64_t> tmp;
		int64_t cc = 0;
		for (size_t i = 0; i < n; i++)
			for (size_t j = 0; j < n; j++) {
				auto d = Uniform1N(gen);
				cc += d;
				tmp.push_back(d);
			}

		size_t idx = 0, m = tmp.size();
		for (size_t i = 0; i < n; i++)
			for (size_t j = 0; j < n; j++) {
				a.add(i, j, tmp[idx]);
				b.add(i, j, tmp[m - idx - 1]);
				idx++;
			}

		PRINT("start solver %d %ld %ld\n", cc, a.balance(), b.balance());
		DOT::Solver solver;

		double dist = solver.bipartite(a, b);

		// TODO: FIX double dist2 = solver.tripartite(a, b);
	}

	return EXIT_SUCCESS;
}
