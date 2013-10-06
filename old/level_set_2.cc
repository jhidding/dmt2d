#include "level_set.h"
#include "conan/mvector.h"
#include "cvector/cvector.h"
#include "conan/boxconfig.h"

#include <set>
#include <vector>
#include <functional>

using namespace Conan;

template <int R>
inline mVector<double, R> c2md(cVector const &x)
{
	mVector<double, R> X;
	for (unsigned i = 0; i < R; ++i) X[i] = x[i];
	return X;
}

Level_set_2::Level_set_2(double const *f)
{
	int const	bits 	= BoxConfig::bits();
	size_t const	size 	= BoxConfig::size();
	cVector dx[2] = { cVector(2, 1, 1), cVector(2, 1, 2) };

	std::set<size_t> C[2];
	for (cVector x(2, bits, 0); x < size; ++x)
	{
		if (f[x] * f[x + dx[0]] < 0.0) C[0].insert(x);
		if (f[x] * f[x + dx[1]] < 0.0) C[1].insert(x);
	}

	while (not C[0].empty())
	{
		cVector x(2, bits, *C[0].begin());
		int direction = 1, i = 0, j = 1;
		Contour_2 L;

		while (true)
		{
			C[i].erase(x);
			cVector target = (direction == 1 ? x : x - dx[j]);
			double frac = - f[x] / (f[x + dx[i]] - f[x]);
			Point p = c2md<2>(x) + c2md<2>(dx[i]) * frac;
			L.push_back(p);

			// quarter turn
			if (C[j].count(target + dx[i]) > 0)
			{
				std::swap(i, j);
				direction = 1;
				x = target + dx[j];
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
			cVector y = (direction == 1 ? x + dx[j] : x - dx[j]);
			if (C[i].count(y) > 0)
			{
				x = y;
				continue;
			}

			break;
		}

		push_back(L);
	}
}

void Level_set_2::output_for_gnuplot(std::ostream &out) const 
{
	std::for_each(begin(), end(),
		[&out] (Contour_2 const &C)
	{
		C.output_for_gnuplot(out);

		out << "\n\n";
	});
}


inline double det2(mVector<double, 2> const &a, mVector<double, 2> const &b)
{
	return a[0] * b[1] - b[0] * a[1];
}

inline double tet2(mVector<double, 2> const &a, mVector<double, 2> const &b, mVector<double, 2> const &c)
{
	return det2(b - a, c - a);
}

Double_level_set_2::Double_level_set_2(double const *f, double const *g)
{
	int const	bits 	= BoxConfig::bits();
	size_t const	size 	= BoxConfig::size();
	cVector dx[2] = { cVector(2, 1, 1), cVector(2, 1, 2) };
	mVector<double, 2> dX[2];
	dX[0] = c2md<2>(dx[0]); dX[1] = c2md<2>(dx[1]);

	std::set<size_t> C;
	for (cVector x(2, bits, 0); x < size; ++x)
	{
		bool f_crosses = 
			   (f[x] * f[x + dx[0]] < 0)
			or (f[x + dx[1]] * f[(x + dx[1]) + dx[0]] < 0)
			or (f[x] * f[x + dx[1]] < 0)
			or (f[x + dx[0]] * f[(x + dx[0]) + dx[1]] < 0);

		bool g_crosses = 
			   (g[x] * g[x + dx[0]] < 0)
			or (g[x + dx[1]] * g[(x + dx[1]) + dx[0]] < 0)
			or (g[x] * g[x + dx[1]] < 0)
			or (g[x + dx[0]] * g[(x + dx[0]) + dx[1]] < 0);

		if (f_crosses and g_crosses)
			C.insert(x);
	}

	std::for_each(C.begin(), C.end(),
		[&] (size_t const &x_)
	{
		cVector x(2, bits, x_);

		std::vector<Point> fc, gc;

		std::function<bool (double const *, cVector const &, int)> 
			p_crosses = 
			[dx] (double const *q, cVector const &r, int i)
		{
			return q[r] * q[r + dx[i]] < 0;
		};

		std::function<Point (double const *, cVector const &, int)> 
			crossing_point =
			[dx, dX] (double const *q, cVector const &r, int i)
		{
			double fraction = - q[r] / (q[r + dx[i]] - q[r]);
			return dX[i] * fraction;
		};

		if (p_crosses(f, x, 0))
			fc.push_back(crossing_point(f, x, 0));
		if (p_crosses(f, x, 1)) 
			fc.push_back(crossing_point(f, x, 1));
		if (p_crosses(f, x + dx[1], 0)) 
			fc.push_back(dX[1] + crossing_point(f, x + dx[1], 0));
		if (p_crosses(f, x + dx[0], 1)) 
			fc.push_back(dX[0] + crossing_point(f, x + dx[0], 1));
		if (p_crosses(g, x, 0)) 	
			gc.push_back(crossing_point(g, x, 0));
		if (p_crosses(g, x, 1)) 	
			gc.push_back(crossing_point(g, x, 1));
		if (p_crosses(g, x + dx[1], 0)) 
			gc.push_back(dX[1] + crossing_point(g, x + dx[1], 0));
		if (p_crosses(g, x + dx[0], 1)) 
			gc.push_back(dX[0] + crossing_point(g, x + dx[0], 1));

		if (gc.size() + fc.size() != 4)
		{
			std::cerr << "#!@#$ \n";
		}
		else if (tet2(fc[0], fc[1], gc[1]) * tet2(fc[0], fc[1], gc[0]) < 0.0)
		{
			Point X = 
				c2md<2>(x) + gc[0] + (gc[1] - gc[0]) *
				det2(fc[1] - fc[0], fc[0] - gc[0]) / 
				det2(fc[1] - fc[0], gc[1] - gc[0]);

			push_back(X);
		}
	});
}

void Double_level_set_2::output_for_gnuplot(double const *ev, std::ostream &out) const 
{
	std::for_each(begin(), end(),
		[&out, ev] (Point const &C)
	{
		size_t i = Conan::BoxConfig::M(C);
		if (ev[i] > 0.0)
			out << C << " " << ev[i] << "\n";
	});
}

void Double_level_set_2::output_for_gnuplot(std::ostream &out) const 
{
	std::for_each(begin(), end(),
		[&out] (Point const &C)
	{
		out << C << "\n";
	});
}

