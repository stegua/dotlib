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

		void resize(size_t t) {
			X.resize(t);
			Y.resize(t);
			B.resize(t);
		}

		void add(int x, int y, double b = 0.0) {
			X.push_back(x);
			Y.push_back(y);
			B.push_back(b);
		}

		bool empty() const {
			return X.empty();
		}

		int getX(size_t i) const { return X[i]; }
		int getY(size_t i) const { return Y[i]; }

		size_t size(void) const { return X.size(); }

		void dump() const {
			for (size_t i = 0, i_max = X.size(); i < i_max; ++i)
				fprintf(stdout, "(%d, %d)\n");
			fprintf(stdout, "\n");
		}

	private:
		// Point coordinates (integers)
		std::vector<int> X;
		std::vector<int> Y;

		// Node balance
		std::vector<double> B;
	};


	class ConvexHull {

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
			int ymax = Ps.getY(0);
			for (size_t i = 0, i_max = Ps.size(); i < i_max; ++i)
				ymax = std::min<int>(ymax, Ps.getY(i));

			std::unordered_map<int, std::vector<int>> Xs;
			for (int i = 0; i < ymax + 1; ++i)
				Xs[Ps.getY(i)].push_back(Ps.getX(i));

			PointCloud2D Bs;
			for (auto& k : Xs) {
				auto vet = k.second;
				std::sort(vet.begin(), vet.end());
				Bs.add(vet.front(), k.first);
				if (vet.size() > 1)
					Bs.add(vet.back(), k.first);
			}


			// Then, filter according to the second axis
			int xmax = Bs.getX(0);
			for (size_t i = 0, i_max = Bs.size(); i < i_max; ++i)
				xmax = std::min<int>(xmax, Bs.getX(i));

			std::unordered_map<int, std::vector<int>> Ys;
			for (int i = 0; i < xmax + 1; ++i)
				Xs[Ps.getX(i)].push_back(Ps.getY(i));


			PointCloud2D Rs;
			for (auto& k : Ys) {
				auto vet = k.second;
				std::sort(vet.begin(), vet.end());
				Bs.add(k.first, vet.front());
				if (vet.size() > 1)
					Bs.add(k.first, vet.back());
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

			PointCloud2D Ss = PolarQuickSort(Cs);
			//		del Ss[Ss.index(anchor)];

			PointCloud2D Hull;
			Hull.add(anchor_x, anchor_y);
			Hull.add(Ss.getX(0), Ss.getY(0));
			size_t cur = 2;
			for (size_t i = 1, i_max = Ss.size(); i < i_max; ++i) {
				while (Det(Hull.getX(cur - 2), Hull.getY(cur - 2),
					Hull.getX(cur - 1), Hull.getY(cur - 1),
					Ss.getX(i), Ss.getY(i)) <= 0) {
					Hull.resize(cur - 1);
					cur--;
					if (Hull.size() < 2)
						break;
				}
				Hull.add(Ss.getX(i), Ss.getY(i));
				cur++;
			}

			return Hull;
		}

	private:
		int anchor_x, anchor_y;
	};

}
