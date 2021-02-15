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
	int n = 512;

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

		a.parse("C:\\Users\\Gualandi\\Google Drive\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data32_1001.csv");
		b.parse("C:\\Users\\Gualandi\\Google Drive\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\ClassicImages\\data32_1002.csv");

		PRINT("start solver %ld %ld\n", a.balance(), b.balance());
		DOT::Solver solver;

		//double dist = solver.solve(a, b);
		//PRINT("ColG => %d: %.6f %.3f sec\n", n, dist, solver.runtime());

		double dist2 = solver.phaseOne(a, b);
		PRINT("ColG => %.6f %.3f sec\n", dist2, solver.runtime());
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

		double dist = solver.solve(a, b);
		PRINT("ColG => %d: %.6f %.3f sec\n", n, dist, solver.runtime());

		double dist2 = solver.phaseOne(a, b);
		PRINT("ColG => %d: %.6f %.3f sec\n", n, dist2, solver.runtime());
	}

	return EXIT_SUCCESS;
}
