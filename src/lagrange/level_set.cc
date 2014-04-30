#include "../base/array.hh"
#include "../base/mvector.hh"
#include "../base/cvector.hh"
#include "../base/boxconfig.hh"
#include <set>
#include <algorithm>

using namespace System;

typedef mVector<double,2> Point;
typedef Array<Point> Contour;
typedef Array<Contour> Level_set;

Level_set compute_level_set_2(ptr<BoxConfig<2>> box, double a, Array<double> f)
{
	Level_set result(0, 0);

	size_t size = box->size();
	cVector<2> ibox(box->bits());
	std::vector<size_t> &dx = ibox.dx_i;

	// each cell has its top and left edges associated with its index.
	// if the levelset crosses an edge add the edge to the collection C[0|1].
	std::set<size_t> C[2];
	for (size_t x; x < size; ++x)
	{
		if ((f[x] - a) * (f[ibox.add(x, dx[0])] - a) < 0.0) C[0].insert(x);
		if ((f[x] - a) * (f[ibox.add(x, dx[1])] - a) < 0.0) C[1].insert(x);
	}

	// while edges crossings are not associated with some contour, continue running
	while (not C[0].empty())
	{
		size_t x = *C[0].begin();
		int direction = 1, i = 0, j = 1;
		Contour L(0);

		while (true)
		{
			// erase last edge from C
			C[i].erase(x);
			// go to the adjacent cell
			size_t target = (direction == 1 ? x : ibox.sub(x, dx[j]));
			// compute the point of intersection, and push it to the levelset
			double frac = - (f[x] - a) / (f[ibox.add(x, dx[i])] - f[x]);
			Point p = ibox.dvec(x) + ibox.dvec(dx[i]) * frac;
			L->push_back(p);

			// quarter turn
			if (C[j].count(ibox.add(target, dx[i])) > 0)
			{
				std::swap(i, j);
				direction = 1;
				x = ibox.add(target, dx[j]);
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
			size_t y = (direction == 1 ? ibox.add(x, dx[j]) : ibox.sub(x, dx[j]));
			if (C[i].count(y) > 0)
			{
				x = y;
				continue;
			}

			break;
		}

		result->push_back(L);
	}

	return result;
}

void print_contour_2(ptr<BoxConfig<2>> box, std::ostream &out, Contour C)
{
	for (unsigned a = 0; a < C.size(); ++a)
	{
		unsigned i = System::modulus(int(a), int(C.size())),
			 j = System::modulus(int(a) - 1, int(C.size()));

		if ((C[i] - C[j]).norm() > box->L() / 2)
			out << "\n\n";

		out << C[i] << "\n";
	}
}

void print_level_set_2(ptr<BoxConfig<2>> box, std::ostream &out, Level_set L)
{
	for (auto C : L) 
	{
		print_contour_2(box, out, C);
		out << "\n\n";
	}
}

