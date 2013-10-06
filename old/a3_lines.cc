#include <set>
#include <functional>
#include <iostream>
#include <algorithm>

#include "level_set.h"
#include "cvector/cvector.h"
#include "kdtree.h"

/*
#include "spline/interpol.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
*/
using namespace Conan;
using namespace std;

template <int R>
inline mVector<double, R> c2md(cVector const &x)
{
	mVector<double, R> X;
	for (unsigned i = 0; i < R; ++i) X[i] = x[i];
	return X;
}

/*
double gsl_function_wrapper(double x, void *params)
{
	std::function<double (double)> *f =
		reinterpret_cast<std::function<double (double)> *>(params);

	return (*f)(x);
}

	const gsl_root_fsolver_type *gsl_T =
		gsl_root_fsolver_bisection;
	gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_T);
	gsl_function F;
	F.function = &gsl_function_wrapper;
*/

template <typename T, int R>
class Distance_to: public kdTree::Distance<Conan::mVector<T, R>, R>
{
	typedef Conan::mVector<T, R> Point;
	Point X;

	public:
		Distance_to(Point const &X_):
			X(X_)
		{}

		virtual double operator()(Point const &Y) const
		{
			return Conan::sqr(X - Y);
		};

		virtual double operator()(kdTree::BoundingBox<Point, R> const &B) const
		{
			bool inside = true;
			double d;

			for (unsigned k = 0; k < R; ++k)
			{
				double U, V, S;
				U = B.min_coord(k) - X[k];
				V = X[k] - B.max_coord(k);
				S = U * V;

				// when inside interval, direction 
				// doesn't add to distance
				if (S < 0)
				{
					d += (U < 0 ? U*U : V*V);
					inside = false;
				}
			}

			return (inside ? -1 : d);
		};
};

A3_lines_2::A3_lines_2(double const *ev,
		double const *H11, double const *H12, double const *H22,
		std::vector<Point> &D3_points, bool versa = true, bool versb = true)
{
	int const	bits 	= BoxConfig::bits();
	size_t const	size 	= BoxConfig::size();
	cVector dx[2] = { cVector(2, 1, 1), cVector(2, 1, 2) };
	mVector<double, 2> dX[2];
	dX[0] = c2md<2>(dx[0]); dX[1] = c2md<2>(dx[1]);

	kdTree::Tree<Point, 2> 
		D3_tree(D3_points.begin(), D3_points.end(),
		[] (Point const &X, int i) -> double { return X[i]; });

	std::function <double (cVector const &, int)>
		grad_ev = [&] (cVector const &x, int i)
			{ return (ev[x + dx[i]] - ev[x - dx[i]]) / 2.0; };

	std::function <double (cVector const &)>
		calc_ia = [&] (cVector const &x)
			{ return grad_ev(x, 0) * (ev[x] - H22[x]) + 
				 grad_ev(x, 1) * (H12[x]); },
		calc_ib = [&] (cVector const &x)
			{ return grad_ev(x, 0) * (H12[x]) + 
				 grad_ev(x, 1) * (ev[x] - H11[x]); };

	std::function<bool (cVector const &, int)> 
		flip_a = [&] (cVector const &x, int i)
			{ return calc_ia(x) * calc_ia(x+ dx[i]) < 0; },
		flip_b = [&] (cVector const &x, int i)
			{ return calc_ib(x) * calc_ib(x + dx[i]) < 0; },
		flip_ea = [&] (cVector const &x, int i)
			{ return ((ev[x] - H22[x]) * (ev[x + dx[i]] - H22[x + dx[i]]) + 
				   H12[x] * H12[x + dx[i]]) < 0; },
		flip_eb = [&] (cVector const &x, int i)
			{ return ((ev[x] - H11[x]) * (ev[x + dx[i]] - H11[x + dx[i]]) + 
				   H12[x] * H12[x + dx[i]]) < 0; };

	std::cerr << "Finding contour candidate points ... \n";
	std::set<size_t> C[2];
	for (cVector x(2, bits, 0); x < size; ++x)
	{
		if (versa)
		{
			if (flip_a(x, 0) and not flip_ea(x, 0))
				C[0].insert(x);

			if (flip_a(x, 1) and not flip_ea(x, 1))
				C[1].insert(x);
		}

		if (versb)
		{
			if (flip_b(x, 0) and not flip_eb(x, 0))
				C[0].insert(x);

			if (flip_b(x, 1) and not flip_eb(x, 1))
				C[1].insert(x);
		}

		/*
		if ((flip_a(x, 0) and not flip_ea(x, 0))
		 or (flip_b(x, 0) and not flip_eb(x, 0)))
			C[0].insert(x);

		if ((flip_a(x, 1) and not flip_ea(x, 1))
		 or (flip_b(x, 1) and not flip_eb(x, 1)))
			C[1].insert(x);*/
	}

	std::cerr << "Tracing contours ... \n";
	while (not C[0].empty())
	{
		cVector q(2, bits, *C[0].begin()), x = q;
		int direction = 1, i = 0, j = 1;
		bool retract = false;
		Line L;

		while (true)
		{
			Point X = c2md<2>(x);
			if (C[i].count(x) > 0) C[i].erase(x);
			cVector target = (direction == 1 ? x : x - dx[j]);

			/* find intersection point using gsl_root solver
			 * on the interpolated function, maybe some overkill.
			 * * * > ------------------------------------------*/
			double fraction = (flip_a(x, i) ?
				-calc_ia(x) / (calc_ia(x + dx[i]) - calc_ia(x)) :
				-calc_ib(x) / (calc_ib(x + dx[i]) - calc_ib(x)));
			L.push_back(X + dX[i] * fraction);

			/* see if we are near a critical point
			 * * * > --------------------------------------------
			Distance_to<double, 2>
				d_to_back(L.back());

			Point Q = D3_tree.nearest_neighbour(d_to_back);
			if (d_to_back(Q) < 1.0)
			{
				L.set_to_have_ends();
				if (retract)
					L.set_end(Q);
				else
					L.set_start(Q);
				L.push_back(Q);
			}
			else
			{*/
				/* decide where to move next, follow contour algorithm
				 * * * > --------------------------------------------*/
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
			//}

			if (not retract)
			{
				direction = -1; i = 0; j = 1; retract = true;
				x = q;
				std::reverse(L.begin(), L.end());
			}
			else 
			{
				/*if (L.unset())
				{*/
					double d = Conan::sqr(L.front() - L.back());
					if (d < 2.0)
						L.set_to_being_loop();
				/*}
				if (L.unset())
				{
					L.set_to_have_ends();
					Distance_to<double, 2>
						d_to_back(L.back()),
						d_to_front(L.front());

					Point Q = D3_tree.nearest_neighbour(d_to_back),
					      U = D3_tree.nearest_neighbour(d_to_front);
					L.set_start(U);
					L.set_end(Q);
				}*/
				break;
			}
		}

		data.push_back(L);
	}

	std::cerr << "\t[done]\n";
}

void A3_lines_2::output_for_gnuplot(double const *ev, std::ostream &out) const
{
	std::for_each(data.begin(), data.end(),
		[&out, ev] (Line const &C)
	{
		C.output_for_gnuplot(ev, out);
		out << "\n\n";
	});
}

#include "conan/interpol.h"

void Critical_line_2::output_for_gnuplot(double const *Eval, std::ostream &out) const
{
	Interpol<double, 2> ev(Eval);

	size_t s = size();
	if (is_loop()) s++;

	bool pen_on_paper = true;
	for (unsigned a = 0; a < s; ++a)
	{
		unsigned i = Conan::modulus(int(a), int(size())),
			 j = Conan::modulus(int(a) - 1, int(size()));

		if (Conan::norm((*this)[i] - (*this)[j]) 
				> Conan::BoxConfig::N() / 2)
			out << "\n\n";

		double v = ev((*this)[i]);
		if (v > 0.0)
		{
			if (not pen_on_paper) out << "\n\n";
			out << (*this)[i] << " " << ev((*this)[i]) << "\n";
			pen_on_paper = true;
		}
		else
		{
			pen_on_paper = false;
		}
	}
}

std::vector<mVector<double, 2>> A3_lines_2::A4_points(double const *Eval)
{
	constexpr int pers = 3;
	std::vector<Point> result;

	Interpol<double, 2> ev(Eval);
	std::for_each(data.begin(), data.end(),
		[&ev, &result] (Line const &L)
	{
		int M = L.size();
		std::vector<double> v;
		std::transform(L.begin(), L.end(), std::back_inserter(v), ev);

		for (int i = pers+2; i < M-pers-2; ++i)
		{
			double sign[2*pers]; bool is_extremum = true;
			for (int j = 0; j < 2*pers; ++j)
			{
				int a0 = Conan::modulus(i + j - pers, M),
				    a1 = Conan::modulus(i + j - pers + 1, M);
				sign[j] = v[a1] - v[a0];
			}
			for (int j = 0; j < pers; ++j)
			{
				if ((sign[j] * sign[2*pers - 1 -j]) > 0)
				{
					is_extremum = false;
					break;
				}
			}

			if (is_extremum)
			{
				result.push_back(L[i]);
			}
		}	
	});

	return result;
}

