/*
 * @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
 *               via Ferrata, 5, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#include <random>

#include "DOT_Solver.h"

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
		std::string SEP = "\\";

		std::string base = "C:\\Users\\Gualandi\\Google "
			"Drive\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";

		DOT::Histogram2D a;
		DOT::Histogram2D b;

		std::string f1 = "data32_1006.csv";
		std::string f2 = "data32_1007.csv";

		a.parse(base + "Shapes" + SEP + f1);
		b.parse(base + "Shapes" + SEP + f2);

		// a.parse("C:"
		//        "\\Users\\gualandi\\Documents\\GitHub\\dotlib\\nearby\\msvc\\Neraby"
		//        "\\x64\\Release\\test1.csv");
		// b.parse("C:"
		//        "\\Users\\gualandi\\Documents\\GitHub\\dotlib\\nearby\\msvc\\Neraby"
		//        "\\x64\\Release\\test2.csv");

		std::string msg = "ClassicImages";

		DOT::Solver solver;

		//solver.tripartiteColgen(a, b);
		//solver.bipartite(a, b);
		solver.colgen(a, b, 20, msg);
	}
	if (false) {
		std::string SEP = "\\";

		std::string base = "C:\\Users\\Gualandi\\Google "
			"Drive\\Ricerca\\DOTA\\data\\DOTmark_1.0\\Data\\";

		std::vector<std::string> dirs = {
			 "ClassicImages", //"Shapes"
			 //"WhiteNoise", "CauchyDensity", "GRFmoderate","MicroscopyImages",
			 //"GRFrough", "GRFsmooth", "LogGRF", "LogitGRF",,
		};

		std::vector<std::string> Fs = {
			 "1001.csv", "1002.csv", "1003.csv", "1004.csv", "1005.csv",
			 "1006.csv", "1007.csv", "1008.csv", "1009.csv", "1010.csv" };

		std::vector<std::string> Ss = { "32", "64", "128", "256", "512" };

		for (const auto& dtype : dirs) {
			for (const auto& S : Ss) {
				for (const auto& f11 : Fs) {
					for (const auto& f22 : Fs)
						if (f11 < f22) {
							DOT::Histogram2D a;
							DOT::Histogram2D b;

							std::string f1 = "data" + S + "_" + f11;
							std::string f2 = "data" + S + "_" + f22;

							a.parse(base + dtype + SEP + f1);
							b.parse(base + dtype + SEP + f2);

							std::string msg = dtype + " " + f11 + " " + f22;

							DOT::Solver solver;

							// for (int tau : {1,  5,  10, 15, 20, 25, 30, 35, 40, 45,
							//                50, 55, 60, 65, 70, 75, 80, 85, 95, 100})
							solver.colgen(a, b, 15, msg);
						}
				}
			}
		}
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
