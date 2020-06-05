/**
 * @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
 *               via Ferrata, 1, I-27100, Pavia, Italy
 *
 * @author stefano.gualandi@gmail.com (Stefano Gualandi)
 *
 */

#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace std {
	template <> struct hash<std::pair<int, int>> {
		inline size_t operator()(const std::pair<int, int>& v) const {
			std::hash<int> int_hasher;
			return int_hasher(v.first) ^ int_hasher(v.second);
		}
	};

}
namespace DOT {

	class PointCloud2D {
	public:
		void remove(int i) {
			std::swap(X[i], X.back());
			std::swap(Y[i], Y.back());
			std::swap(B[i], B.back());
			X.resize(X.size() - 1);
			Y.resize(Y.size() - 1);
			B.resize(B.size() - 1);
		}

		void remove(int x, int y) {
			auto p = std::make_pair(x, y);
			if (M.find(p) != M.end()) {
				size_t i = M.at(p);
				std::swap(X[i], X.back());
				std::swap(Y[i], Y.back());
				std::swap(B[i], B.back());
				X.pop_back();
				Y.pop_back();
				B.pop_back();
				M.erase(p);
			}
		}

		void reserve(size_t t) {
			X.reserve(t);
			Y.reserve(t);
			B.reserve(t);
		}

		void pop_back() {
			remove(X.back(), Y.back());
		}

		void shrink_to_fit() {
			X.shrink_to_fit();
			Y.shrink_to_fit();
			B.shrink_to_fit();
		}

		void add(int x, int y, double b = 0.0) {
			auto p = std::make_pair(x, y);
			if (M.find(p) == M.end()) {
				M[p] = X.size();
				X.push_back(x);
				Y.push_back(y);
				B.push_back(b);
			}
		}

		void setX(size_t i, int x) {
			X[i] = x;
		}
		void setY(size_t i, int y) {
			Y[i] = y;
		}
		void setB(size_t i, double b) {
			B[i] = b;
		}

		// Merge all the points contained in "other" into this object.
		// The node balance are taken from the "other" object.
		void merge(const PointCloud2D& other) {
			std::unordered_map< std::pair<int, int>, size_t > O = other.getM();

			for (const auto& p : O) {
				size_t j = p.second;
				if (M.find(p.first) != M.end()) {
					size_t i = M.at(p.first);
					B[i] = other.getB(j);
				}
				else {
					fprintf(stdout, "point missing: %d %d\n", p.first.first, p.first.second);
				}
			}

		}

		void append(const PointCloud2D& other) {
			for (size_t i = 0, i_max = other.size(); i < i_max; ++i)
				add(other.getX(i), other.getY(i), other.getB(i));
		}

		double balance() {
			double t = 0;
			for (size_t i = 0, i_max = B.size(); i < i_max; ++i)
				t += B[i];
			return t;
		}


		bool empty() const {
			return X.empty();
		}

		int getX(size_t i) const { return X[i]; }
		int getY(size_t i) const { return Y[i]; }
		double getB(size_t i) const { return B[i]; }

		size_t size(void) const { return X.size(); }

		void dump(const std::string& msg = "") const {
			if (!msg.empty())
				fprintf(stdout, "%s\n", msg.c_str());
			for (size_t i = 0, i_max = X.size(); i < i_max; ++i)
				fprintf(stdout, "(%d, %d, %f)\n", X[i], Y[i], B[i]);
			fprintf(stdout, "\n");
			fflush(stdout);
		}

		const std::unordered_map< std::pair<int, int>, size_t >& getM() const { return M; }

	private:
		// Point coordinates (integers)
		std::vector<int> X;
		std::vector<int> Y;
		// Pair to indices
		std::unordered_map< std::pair<int, int>, size_t > M;
		// Node balance
		std::vector<double> B;

	};

	// Parse data from file, with format: i j b1 b1
	PointCloud2D parse(const std::string& filename, char sep = ' ', int off = 0) {
		std::ifstream in_file(filename);

		if (!in_file) {
			fprintf(stderr, "FATAL ERROR: Cannot open file %s", filename.c_str());
			exit(EXIT_FAILURE);
		}

		PointCloud2D Rs;
		std::vector<double> Bs;
		std::string         line;

		// Read first line
		double tot_a = 0;
		double tot_b = 0;
		while (std::getline(in_file, line)) {
			std::stringstream   lineStream(line);
			std::string         cell;

			std::getline(lineStream, cell, sep);
			int x = std::stoi(cell);
			std::getline(lineStream, cell, sep);
			int y = std::stoi(cell);
			std::getline(lineStream, cell, sep);
			double a = std::stof(cell);
			std::getline(lineStream, cell, sep);
			double b = std::stof(cell);

			if (fabs(a - b) > 1e-10) {
				tot_a += a;
				tot_b += b;
				// Check if grid start in position 1 or 0 with parameter "off"
				Rs.add(x - off, y - off, a);
				Bs.emplace_back(b);
			}
		};

		// Release resource as soon as possible
		in_file.close();
		// Use as few memory as possible
		Rs.shrink_to_fit();
		// normalize data (rescaling)
		if (Rs.size() != Bs.size())
			exit(0);
		double tot = 0;
		for (size_t i = 0, i_max = Bs.size(); i < i_max; ++i) {
			tot += Rs.getB(i) / tot_a - Bs[i] / tot_b;
			Rs.setB(i, Rs.getB(i) / tot_a - Bs[i] / tot_b);
		}

		fprintf(stdout, "tot flow: %f\n", tot);

		return Rs;
	}

	// Class for computing the convex hull
	class ConvexHull {
	public:
		// Compute polar among between two points
		double PolarAngle(int ax, int ay, int bx = -1, int by = -1) const {
			int cx = bx, cy = by;
			if (bx == -1) {
				cx = anchor_x;
				cy = anchor_y;
			}
			int x_span = ax - cx;
			int y_span = ay - cy;
			return atan2(y_span, x_span);
		}

		// Square Euclidean distance
		int Distance(int ax, int ay, int bx = -1, int by = -1) const {
			int cx = bx, cy = by;
			if (bx == -1) {
				cx = anchor_x;
				cy = anchor_y;
			}
			int x_span = ax - cx;
			int y_span = ay - cy;
			return pow(y_span, 2) + pow(x_span, 2);
		}

		// Determinant to detect direction
		int Det(int ax, int ay, int bx, int by, int cx, int cy) const {
			return  (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
		}


		PointCloud2D PolarQuickSort(PointCloud2D& Ls) {
			if (Ls.size() <= 1)
				return Ls;

			PointCloud2D smaller, equal, larger;

			double pivot_ang = PolarAngle(Ls.getX(0), Ls.getY(0));

			for (size_t i = 0, i_max = Ls.size(); i < i_max; ++i) {
				double p_ang = PolarAngle(Ls.getX(i), Ls.getY(i));
				if (p_ang < pivot_ang) {
					smaller.add(Ls.getX(i), Ls.getY(i));
				}
				else {
					if (p_ang == pivot_ang)
						equal.add(Ls.getX(i), Ls.getY(i));
					else
						larger.add(Ls.getX(i), Ls.getY(i));
				}
			}

			auto l1 = PolarQuickSort(smaller);
			while (!equal.empty()) {
				size_t min_idx = 0;
				for (size_t i = 0, i_max = equal.size(); i < i_max; ++i)
					if (Distance(equal.getX(i), equal.getY(i)))
						min_idx = i;
				l1.add(equal.getX(min_idx), equal.getY(min_idx));
				equal.remove(min_idx);
			}
			auto l3 = PolarQuickSort(larger);
			for (size_t i = 0, i_max = l3.size(); i < i_max; ++i)
				l1.add(l3.getX(i), l3.getY(i));
			return l1;
		}


		// Filter the main point along a given direction
		PointCloud2D FilterAxis(const PointCloud2D& Ps) {
			// First filter along an axis
			std::unordered_map<int, std::vector<int>> Xs;
			for (size_t i = 0, i_max = Ps.size(); i < i_max; ++i) {
				int key = Ps.getY(i);
				if (Xs.find(key) == Xs.end()) {
					std::vector<int> tmp;
					Xs[key] = tmp;
				}
				Xs.at(key).push_back(Ps.getX(i));
			}

			PointCloud2D Bs;
			for (auto& k : Xs) {
				auto vet = k.second;
				std::sort(vet.begin(), vet.end());
				Bs.add(vet.front(), k.first);
				if (vet.size() > 1)
					Bs.add(vet.back(), k.first);
			}

			// Then, filter according to the second axis
			std::unordered_map<int, std::vector<int>> Ys;
			for (size_t i = 0, i_max = Bs.size(); i < i_max; ++i) {
				int key = Bs.getX(i);
				if (Ys.find(key) == Ys.end()) {
					std::vector<int> tmp;
					Ys[key] = tmp;
				}
				Ys.at(key).push_back(Bs.getY(i));
			}

			PointCloud2D Rs;
			for (auto& k : Ys) {
				auto vet = k.second;
				std::sort(vet.begin(), vet.end());
				Rs.add(k.first, vet.front());
				if (vet.size() > 1)
					Rs.add(k.first, vet.back());
			}

			return Rs;
		}

		// Find convex hull of given set of points
		PointCloud2D find(const PointCloud2D& Ps) {
			// Preprocessing
			auto Cs = FilterAxis(Ps);
			// Find anchor point
			int min_idx = -1;
			for (size_t i = 0, i_max = Cs.size(); i < i_max; ++i) {
				if (min_idx == -1 || Cs.getY(i) < Cs.getY(min_idx))
					min_idx = i;
				if (Cs.getY(i) == Cs.getY(min_idx) && Cs.getX(i) < Cs.getX(min_idx))
					min_idx = i;
			}

			anchor_x = Cs.getX(min_idx);
			anchor_y = Cs.getY(min_idx);
			Cs.remove(anchor_x, anchor_y);

			PointCloud2D Ss = PolarQuickSort(Cs);

			PointCloud2D Hull;
			Hull.add(anchor_x, anchor_y);
			Hull.add(Ss.getX(0), Ss.getY(0));
			for (size_t i = 1, i_max = Ss.size(); i < i_max; ++i) {
				int cur = Hull.size();
				while (Det(Hull.getX(cur - 2), Hull.getY(cur - 2),
					Hull.getX(cur - 1), Hull.getY(cur - 1),
					Ss.getX(i), Ss.getY(i)) <= 0) {
					Hull.pop_back();
					if (Hull.size() < 2)
						break;
					cur = Hull.size();
				}
				Hull.add(Ss.getX(i), Ss.getY(i));
			}

			// Add all points along the convex hull
			PointCloud2D Rs;
			size_t n = Hull.size();
			for (size_t i = 0, i_max = n - 1; i < i_max; i++) {
				auto tmp = WalkGrid(Hull.getX(i), Hull.getY(i), Hull.getX(i + 1), Hull.getY(i + 1));
				Rs.append(tmp);
			}
			auto tmp = WalkGrid(Hull.getX(n - 1), Hull.getY(n - 1), Hull.getX(0), Hull.getY(0));
			Rs.append(tmp);
			return Rs;
		}

		// Find all the points connecting two dots
		PointCloud2D WalkGrid(int ax, int ay, int bx, int by) {
			int dx = bx - ax;
			int dy = by - ay;

			if (dx == 0) {
				PointCloud2D ps;
				for (int i = std::min(ay, by), i_max = std::max(ay, by); i < i_max; ++i)
					ps.add(ax, i);
				return ps;
			}

			if (dy == 0) {
				PointCloud2D ps;
				for (int i = std::min(ax, bx), i_max = std::max(ax, bx); i < i_max; ++i)
					ps.add(i, ay);
				return ps;
			}

			int nx = abs(dx);
			int ny = abs(dy);

			int sign_x = (dx > 0 ? 1 : -1);
			int sign_y = (dy > 0 ? 1 : -1);

			int px = ax;
			int py = ay;
			PointCloud2D points;
			points.add(px, py);

			int ix = 0, iy = 0;
			while (ix < nx || iy < ny) {
				if ((0.5 + ix) / nx < (0.5 + iy) / ny) {
					px += sign_x;
					ix += 1;
				}
				else {
					py += sign_y;
					iy += 1;
				}
				points.add(px, py);
			}

			return points;
		}

		// Find all the point in the interior of the convex hull
		PointCloud2D FillHull(const PointCloud2D& Ps) const {
			int x_max = -1;
			for (size_t i = 0, i_max = Ps.size(); i < i_max; ++i)
				x_max = std::max(x_max, Ps.getX(i));

			std::vector<int> Xmin(x_max + 1, std::numeric_limits<int>::max());
			std::vector<int> Xmax(x_max + 1, -1);

			for (size_t i = 0, i_max = Ps.size(); i < i_max; ++i) {
				int x = Ps.getX(i);
				Xmin[x] = std::min(Xmin[x], Ps.getY(i));
				Xmax[x] = std::max(Xmax[x], Ps.getY(i));
			}

			for (int i = 0; i < x_max + 1; ++i)
				if (Xmin[i] == std::numeric_limits<int>::max() || Xmax[i] == -1)
					exit(-1);

			PointCloud2D Rs;
			for (int x = 0; x < x_max + 1; ++x)
				for (size_t y = Xmin[x], y_max = Xmax[x] + 1; y < y_max; ++y)
					Rs.add(x, y);

			return Rs;
		}

	private:
		int anchor_x, anchor_y;
	};

}
