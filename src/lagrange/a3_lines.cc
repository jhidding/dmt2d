#include "../base/mvector.hh"
#include "../base/cvector.hh"
#include "../base/array.hh"
#include "../base/misc.hh"
#include "../base/boxconfig.hh"

#include "../misc/hessian.hh"

#include <functional>
#include <set>
#include <iostream>

//#include "../misc/kdtree.h"
//#include "../misc/interpol.h"

using namespace System;

typedef mVector<double,2> Point;
typedef Array<Point> Line;
typedef Array<Line> A3_lines;

A3_lines compute_a3_lines_2(ptr<BoxConfig<2>> box, Array<double> ev, DMT::Hessian<2> H, Array<Point> D4)
{
	// setup some variables
	A3_lines result(0, 0);
	cVector<2> b(box->bits());
	size_t size = box->size();
	std::vector<size_t> &dx = b.dx_i;
	Point dX[2]; for (unsigned i = 0; i < 2; ++i) { dX[i] = b.dvec(dx[i]); }

	// make a kdTree of the D4 points, so we can search efficiently
	//kdTree::Tree<Point, 2>
	//	D4_tree(D4.begin(), D4.end(),
	//	[] (Point const &X, int i) -> double { return X[i]; });

	// gradient of the eigenvalue
	std::function <double (size_t, int)>
		grad_ev = [&] (size_t x, int i)
			{ return (ev[b.add(x, dx[i])] - ev[b.sub(x, dx[i])]) / 2.0; };

	// there are two functions to calculate the direction of the
	// eigenvector if one is in local conflict due to a topological defect,
	// we can use the other.
	std::function <double (size_t)>
		calc_ia = [&] (size_t x)
			{ return grad_ev(x, 0) * (ev[x] - H[2][x]) + 
				 grad_ev(x, 1) * (H[1][x]); },
		calc_ib = [&] (size_t x)
			{ return grad_ev(x, 0) * (H[1][x]) + 
				 grad_ev(x, 1) * (ev[x] - H[0][x]); };

	std::function<bool (size_t, int)> 
		flip_a = [&] (size_t x, int i)
			{ return calc_ia(x) * calc_ia(b.add(x, dx[i])) < 0; },
		flip_b = [&] (size_t x, int i)
			{ return calc_ib(x) * calc_ib(b.add(x, dx[i])) < 0; },
		flip_ea = [&] (size_t x, int i)
			{ return ((ev[x] - H[2][x]) * (ev[b.add(x, dx[i])] - H[2][b.add(x, dx[i])]) + 
				   H[1][x] * H[1][b.add(x, dx[i])]) < 0; },
		flip_eb = [&] (size_t x, int i)
			{ return ((ev[x] - H[0][x]) * (ev[b.add(x, dx[i])] - H[0][b.add(x, dx[i])]) + 
				   H[1][x] * H[1][b.add(x, dx[i])]) < 0; };

	std::cerr << "Finding contour candidate points ... \n";
	std::set<size_t> C[2];
	for (size_t x = 0; x < size; ++x)
	{
		if (flip_a(x, 0) and not flip_ea(x, 0))
			C[0].insert(x);

		if (flip_a(x, 1) and not flip_ea(x, 1))
			C[1].insert(x);
	}

	std::cerr << "Tracing contours ... \n";
	while (not C[0].empty())
	{
		size_t q = *C[0].begin(), x = q;
		int direction = 1, i = 0, j = 1;
		bool retract = false;
		Line L(0);

		while (true)
		{
			Point X = b.dvec(x);
			if (C[i].count(x) > 0) C[i].erase(x);
			size_t target = (direction == 1 ? x : b.sub(x, dx[j]));

			double fraction = (flip_a(x, i) ?
				-calc_ia(x) / (calc_ia(x + dx[i]) - calc_ia(x)) :
				-calc_ib(x) / (calc_ib(x + dx[i]) - calc_ib(x)));
			L->push_back(X + dX[i] * fraction);

			// quarter turn
			if (C[j].count(b.add(target, dx[i])) > 0)
			{
				std::swap(i, j);
				direction = 1;
				x = b.add(target, dx[j]);
				continue;
			}

			// quarter turn other way
			if (C[j].count(target) > 0)
			{
				std::swap(i, j);
				direction = -1;
				x = target;
				continue;
			}

			// ahead
			size_t y = (direction == 1 ? b.add(x, dx[j]) : b.sub(x, dx[j]));
			if (C[i].count(y) > 0)
			{
				x = y;
				continue;
			}

			if (not retract)
			{
				direction = -1; i = 0; j = 1; retract = true;
				x = q;
				std::reverse(L.begin(), L.end());
			}
			else 
			{
				break;
			}
		}

		result->push_back(L);
	}

	std::cerr << "\t[done]\n";
	return result;
}

