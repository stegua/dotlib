/*
* @fileoverview Copyright (c) 2017-20, Stefano Gualandi,
*               via Ferrata, 1, I-27100, Pavia, Italy
*
* @author stefano.gualandi@gmail.com (Stefano Gualandi)
*
*/

// FROM STDLIB
#include <cassert>
#include <cmath>

#include <string>
#include <vector>

#include <fstream>
#include <sstream>
#include <iostream>

#include <algorithm>
#include <numeric>
#include <chrono>

// In order to use PRId64 and PRIu64
#include <inttypes.h>

// Args parser library
#include <args.hxx>

// Yocta Logger library
#include <yocta_logger.hh>

// My Network Simplex
#include "DOT_NetSimplex.h"
#include "DOT_ConvexHull.h"


// Global variable logger
namespace DOT {
	yocta::Logger logger;
}

using namespace DOT;


/**
* @brief Two dimensional histogram
*/
class Histogram2D {
public:
	// Std c'tor
	Histogram2D(const std::string& filename) {
		readFromFile(filename);
	}

	// Rule of five?
	// ...

	// Parse from file
	void readFromFile(const std::string& filename) {
		std::ifstream in_file(filename);

		if (!in_file) {
			logger.error("FATAL ERROR: Cannot open file %s", filename.c_str());
			exit(EXIT_FAILURE);
		}

		// Read first line
		auto read_row = [&](size_t i) {
			size_t j = 0;
			std::string         line;
			std::getline(in_file, line);
			std::stringstream   lineStream(line);
			std::string         cell;

			while (std::getline(lineStream, cell, ',')) {
				data.push_back(stod(cell));
				++j;
			}

			return j;
		};

		// Read first row, and return row length
		n = read_row(0);

		for (size_t i = 1; i < n; ++i)
			read_row(i);

		// Use as few memory as possible
		data.shrink_to_fit();

		in_file.close();
	}

	// Get an elelment. Index data as a matrix
	double get(size_t x, size_t y) const {
		return data[x * n + y];
	};

	// Get dimension of histogram
	size_t getN() const {
		return n;
	}

	// Dump file to stdout
	void dump() const {
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < n; ++j)
				fprintf(stdout, "%f ", get(i, j));
			fprintf(stdout, "\n");
		}
	}

	// Sum over all bins
	double computeTotalWeight() const {
		double t = 0;
		for (double v : data)
			t += v;
		return t;
	}

private:
	size_t n;   // Histogram of size n*n
	std::vector<double>  data;   // Histogram data a single array (contiguos in memory)
};


// In standarc C++17 it is better to use function std::gcd
int gcd(int _a, int _b) {
	int a = (_a >= 0 ? _a : -_a);
	int b = (_b >= 0 ? _b : -_b);
	while (b != 0) {
		int t = b;
		b = a % b;
		a = t;
	}
	return a;
}

#ifdef _WIN32
const auto& GCD = static_cast<int(*)(int, int)>(std::gcd);
#else
#include <algorithm>
const auto& GCD = gcd;/// static_cast<int(*)(int, int)>(std::__gcd);
#endif


// Default test: Wasserstein distance order 1, ground distance L2, two classic images
void solver_W1L2(const std::string& f1, const std::string& f2, int L) {

	logger.note("Parameter L = %d", L);

	// List of pair of coprimes number between (-L, L)
	std::vector<std::pair<int, int>> coprimes;
	for (int v = -L; v <= L; ++v)
		for (int w = -L; w <= L; ++w)
			if (!(v == 0 && w == 0) && GCD(v, w) == 1)
				coprimes.emplace_back(v, w);
	// Use as few memory as possible
	coprimes.shrink_to_fit();

	logger.note("Computed: %d coprimes", (int)coprimes.size());

	// Readin input histograms
	Histogram2D h1(f1);
	Histogram2D h2(f2);

	int s = h1.getN();
	int n = s * s;
	int Co = coprimes.size();

	logger.note("Histogram dimension: s=%d, n=%d", s, n);

	auto ID = [&s](int x, int y) {
		return x * s + y;
	};

	typedef double FlowType;
	typedef double CostType;

	// Build the graph for min cost flow
	NetSimplex<FlowType, CostType> simplex(n, n * Co);

	// add first d source nodes
	for (int i = 0; i < s; ++i)
		for (int j = 0; j < s; ++j)
			simplex.addNode(ID(i, j), h1.get(i, j) - h2.get(i, j));

	int m = 0;
	for (int i = 0; i < s; ++i)
		for (int j = 0; j < s; ++j)
			for (const auto& p : coprimes) {
				int v = p.first;
				int w = p.second;
				if (i + v >= 0 && i + v < s && j + w >= 0 && j + w < s) {
					simplex.setArc(m++, ID(i, j), ID(i + v, j + w), sqrt(pow(v, 2) + pow(w, 2)));
				}
			}

	logger.note("Network, ready. Size: n=%d, m=%d", n, m);

	// Time vars
	using namespace std;
	using namespace std::chrono;
	// Start time
	time_point<system_clock> start = system_clock::now();


	// Solve the problem to compute the distance
	NetSimplex<FlowType, CostType>::ProblemType status = simplex.run();

	double distance = std::numeric_limits<CostType>::max();
	if (status != NetSimplex<FlowType, CostType>::INFEASIBLE &&
		status != NetSimplex<FlowType, CostType>::UNBOUNDED)
		distance = simplex.totalCost();

	// End time
	time_point<system_clock> end = system_clock::now();
	duration<double> inlineTimeElapsed = end - start;

	logger.note("Distance: %f, Time: %ld ms",
		distance, duration_cast<milliseconds>(inlineTimeElapsed).count());
}

// Default test: Wasserstein distance order 1, ground distance L2, two classic images
void solver_W1L2_hull(const std::string& filename, int L) {
	logger.note("Parameter L = %d", L);

	// List of pair of coprimes number between (-L, L)
	std::vector<std::pair<int, int>> coprimes;
	for (int v = -L; v <= L; ++v)
		for (int w = -L; w <= L; ++w)
			if (!(v == 0 && w == 0) && GCD(v, w) == 1)
				coprimes.emplace_back(v, w);
	// Use as few memory as possible
	coprimes.shrink_to_fit();
	int Co = coprimes.size();

	logger.note("Computed: %d coprimes", (int)coprimes.size());

	// Read input data
	PointCloud2D ps = parse(filename);
	logger.note("File in memory, with %d bins. Start convex hull computation", (int)ps.size());
	// Compute convex hull
	ConvexHull ch;
	PointCloud2D As = ch.find(ps);
	PointCloud2D Rs = ch.FillHull(As);
	Rs.merge(ps);
	int n0 = ps.size();
	int n = Rs.size();

	logger.note("Input measures have %d points, convex hull %d points", n0, n);

	// Compute xmin, xmax, ymin, ymax for each axis
	int xmax = 0;
	int ymax = 0;
	for (int i = 0; i < n; ++i) {
		xmax = std::max(xmax, Rs.getX(i));
		ymax = std::max(ymax, Rs.getY(i));
	}
	xmax++;
	ymax++;

	logger.note("Histogram within coordinates: xmax=%d, ymax=%d", xmax, ymax);

	// Binary vector for positions
	auto ID = [&ymax](int x, int y) {
		return x * ymax + y;
	};

	std::vector<bool> M(xmax * ymax, false);
	for (int i = 0; i < n; ++i)
		M[ID(Rs.getX(i), Rs.getY(i))] = true;

	std::vector<int> H(xmax * ymax, 0);
	for (int i = 0; i < n; ++i)
		H[ID(Rs.getX(i), Rs.getY(i))] = i;

	typedef double FlowType;
	typedef double CostType;

	// Build the graph for min cost flow
	NetSimplex<FlowType, CostType> simplex(n, n * Co);

	// add first d source nodes
	for (int i = 0; i < n; ++i)
		simplex.addNode(i, Rs.getB(i));

	int m = 0;
	for (int h = 0; h < n; ++h) {
		int i = Rs.getX(h);
		int j = Rs.getY(h);
		for (const auto& p : coprimes) {
			int v = p.first;
			int w = p.second;
			if (i + v >= 0 && i + v < xmax && j + w >= 0 && j + w < ymax && M[ID(i + v, j + w)]) {
				int ff = H[ID(i + v, j + w)];
				simplex.setArc(m++, h, ff, sqrt(pow(v, 2) + pow(w, 2)));
			}
		}
	}

	logger.note("Network, ready. Size: n=%d, m=%d", n, m);

	// Time vars
	using namespace std;
	using namespace std::chrono;
	// Start time
	time_point<system_clock> start = system_clock::now();

	// Solve the problem to compute the distance
	NetSimplex<FlowType, CostType>::ProblemType status = simplex.run();

	double distance = std::numeric_limits<CostType>::max();
	if (status != NetSimplex<FlowType, CostType>::INFEASIBLE &&
		status != NetSimplex<FlowType, CostType>::UNBOUNDED)
		distance = simplex.totalCost();

	// End time
	time_point<system_clock> end = system_clock::now();
	duration<double> inlineTimeElapsed = end - start;

	logger.note("Distance: %f, Time: %ld ms",
		distance, duration_cast<milliseconds>(inlineTimeElapsed).count());
}

// Main entry point
#ifdef CONVEXHULL

int main(int argc, char* argv[]) {
	if (false) {
		PointCloud2D ps;
		ps.add(2, 3);
		ps.add(2, 4);
		ps.add(3, 4);
		ps.add(3, 5);
		ps.add(3, 6);
		ps.add(3, 7);
		ps.add(1, 4);
		ps.add(3, 9);
		ps.add(1, 6);

		ConvexHull ch;
		PointCloud2D rs = ch.find(ps);

		rs.dump();
	}
	else {
		PointCloud2D ps = parse("data/convex_test.csv");
		ps.dump("prima");
		ConvexHull ch;
		PointCloud2D rs = ch.find(ps);
		rs.merge(ps);
		rs.dump("dopo");
	}
}

#else
int main(int argc, char* argv[]) {
	args::ArgumentParser parser("DOTLib v0.5.0: Discrete Optimal Transport library,\nby Stefano Gualandi (stefano.gualandi@gmail.com), 2019-2020.", "");
	args::HelpFlag help(parser, "help", "Display this help menu", { 'h', "help" });
	args::ValueFlag<std::string> h1_filename(parser, "h1_filename", "filename of the first histogram", { 'f', "h1" });
	args::ValueFlag<std::string> h2_filename(parser, "h2_filename", "filename of the second histogram", { 'g', "h2" });
	args::ValueFlag<std::string> s_filename(parser, "s_filename", "filename of spatial data", { 's', "s1" });
	args::ValueFlag<std::string> log_filename(parser, "log_filename", "logger filename", { 'l', "log" });

	args::ValueFlag<int> side(parser, "side", "Dimension of the subgrid size L", { 'L', "sides" });
	args::ValueFlag<int> verbosity(parser, "verbosity", "Verbosity level of logger", { 'v', "ver" });

	args::Group nomrs(parser, "This group is all exclusive:", args::Group::Validators::AtMostOne);
	args::Flag norm2(nomrs, "L2", "Norm-2 ground distance", { 'b', "l_2" });

	try {
		if (argc == 1) {
			std::cout << argv[0] << parser;
			return EXIT_SUCCESS;
		}

		parser.ParseCLI(argc, argv);

	}
	catch (args::Help) {
		std::cout << parser;
		return EXIT_SUCCESS;
	}
	catch (args::ParseError e) {
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return EXIT_FAILURE;
	}
	catch (args::ValidationError e) {
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return EXIT_FAILURE;
	}

	// Start the command line options
	try {
		// Set logger first
		if (log_filename)
			logger.setFileStream(args::get(log_filename));
		if (verbosity)
			logger.setVerbosityLevel(static_cast<yocta::VerbosityLevel>(args::get(verbosity)));

		// For W1L2 Parameter L
		int L = 2;
		if (side)
			L = args::get(side);

		if (s_filename) {
			std::string s1 = args::get(s_filename);
			solver_W1L2_hull(s1, L);
		}
		else {
			std::string f1 = "data32_1001.csv";
			std::string f2 = "data32_1002.csv";
			if (h1_filename && h2_filename) {
				f1 = args::get(h1_filename);
				f2 = args::get(h2_filename);
			}
			// Run single test
			solver_W1L2(f1, f2, L);
		}

		exit(EXIT_SUCCESS);
	}
	catch (std::exception e) {
		logger.error("Runtime error: %s", e.what());
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

#endif // TEST_CONVEX_HULL
