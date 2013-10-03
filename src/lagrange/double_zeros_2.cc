#include <functional>
#include <set>
#include <vector>
#include "../base/array.hh"
#include "../base/boxconfig.hh"
#include "../base/mvector.hh"
#include "../base/cvector.hh"

using namespace System;

inline double det2(mVector<double, 2> const &a, mVector<double, 2> const &b)
{
	return a[0] * b[1] - b[0] * a[1];
}

inline double tet2(mVector<double, 2> const &a, mVector<double, 2> const &b, mVector<double, 2> const &c)
{
	return det2(b - a, c - a);
}

Array<mVector<double, 2>> compute_double_zeros_2(ptr<BoxConfig<2>> box, Array<double> f, Array<double> g)
{
	typedef mVector<double,2> Point;
	Array<Point> result(0);

	size_t size = box->size();
	cVector<2> ibox(box->bits());
	std::vector<size_t> &dx = ibox.dx_i;

	Point dX[2];
	dX[0] = ibox.dvec(dx[0]); dX[1] = ibox.dvec(dx[1]);

	// search points where both <f> and <g> have roots, and collect them in <C>
	std::set<size_t> C;
	for (size_t x; x < size; ++x)
	{
		bool f_crosses = 
			   (f[x] * f[ibox.add(x, dx[0])] < 0)
			or (f[ibox.add(x, dx[1])] * f[ibox.add(ibox.add(x, dx[1]), dx[0])] < 0)
			or (f[x] * f[ibox.add(x, dx[1])] < 0)
			or (f[ibox.add(x, dx[0])] * f[ibox.add(ibox.add(x, dx[0]), dx[1])] < 0);

		bool g_crosses = 
			   (g[x] * g[ibox.add(x, dx[0])] < 0)
			or (g[ibox.add(x, dx[1])] * g[ibox.add(ibox.add(x, dx[1]), dx[0])] < 0)
			or (g[x] * g[ibox.add(x, dx[1])] < 0)
			or (g[ibox.add(x, dx[0])] * g[ibox.add(ibox.add(x, dx[0]), dx[1])] < 0);

		if (f_crosses and g_crosses)
			C.insert(x);
	}

	std::function<bool (Array<double>,  size_t, int)> 
		p_crosses = [&ibox] (Array<double> q, size_t r, int i)
	{
		return q[r] * q[ibox.add(r, ibox.dx_i[i])] < 0;
	};

	std::function<Point (Array<double>, size_t, int)> 
		crossing_point = [&ibox, dX] (Array<double> q, size_t r, int i)
	{
		double fraction = - q[r] / (q[ibox.add(r, ibox.dx_i[i])] - q[r]);
		return dX[i] * fraction;
	};

	// loop over all points in <C> and find accurate location of double root,
	// by linear interpolation of both functions.
	for (auto x : C)
	{
		std::vector<Point> fc, gc;

		if (p_crosses(f, x, 0))
			fc.push_back(crossing_point(f, x, 0));
		if (p_crosses(f, x, 1)) 
			fc.push_back(crossing_point(f, x, 1));
		if (p_crosses(f, ibox.add(x, dx[1]), 0)) 
			fc.push_back(dX[1] + crossing_point(f, ibox.add(x, dx[1]), 0));
		if (p_crosses(f, ibox.add(x, dx[0]), 1)) 
			fc.push_back(dX[0] + crossing_point(f, ibox.add(x, dx[0]), 1));
		if (p_crosses(g, x, 0)) 	
			gc.push_back(crossing_point(g, x, 0));
		if (p_crosses(g, x, 1)) 	
			gc.push_back(crossing_point(g, x, 1));
		if (p_crosses(g, ibox.add(x, dx[1]), 0)) 
			gc.push_back(dX[1] + crossing_point(g, ibox.add(x, dx[1]), 0));
		if (p_crosses(g, ibox.add(x, dx[0]), 1)) 
			gc.push_back(dX[0] + crossing_point(g, ibox.add(x, dx[0]), 1));

		if (gc.size() + fc.size() != 4)
		{
			std::cerr << "#!@#$ \n";
		}
		else if (tet2(fc[0], fc[1], gc[1]) * tet2(fc[0], fc[1], gc[0]) < 0.0)
		{
			Point X = 
				ibox.dvec(x) + gc[0] + (gc[1] - gc[0]) *
				det2(fc[1] - fc[0], fc[0] - gc[0]) / 
				det2(fc[1] - fc[0], gc[1] - gc[0]);

			result->push_back(X);
		}
	}
}

