#pragma once
#include "../base/system.hh"
#include <queue>

namespace Misc
{
	namespace Interpol {

	using System::cVector;
	using System::mVector;

	template <typename Q, unsigned R>
	class Linear
	{
		using dVector = System::mVector<double, R>;
		using iVector = System::mVector<int, R>;

		System::ptr<System::BoxConfig<R>> box;
		Q						f;

		public:
			typedef typename Q::value_type FT;
			typedef mVector<double, R> Point;

			Linear(System::ptr<System::BoxConfig<R>> box_, Q f_):
				box(box_), f(f_)
			{}

			typename Q::value_type operator()(Point const &x) const
			{
				iVector origin;
				dVector A[2];
				//double r = box->scale();

				for (unsigned k = 0; k < R; ++k)
				{
					origin[k] = (int)(x[k]);
					A[1][k] = x[k] - origin[k];
					A[0][k] = 1 - A[1][k];
				}

				auto s = box->m2c(origin);
				typename Q::value_type v(0);

				for (auto dx : box->sq())
				{
					double z = 1;
					for (unsigned k = 0; k < R; ++k) 
						z *= A[dx[k]][k];

					v += f[s + dx] * z;
				}

				return v;
			}
	};

	template <typename F, int R>
	class Spline {};

	template <typename F>
	class Spline<F, 1>
	{
		cVector<2>	loop;
		enum { R = 1 };
		cVector<1>	box;

		static double	A[16];
		F		fd;

		typedef mVector<double, R> 	fVector;
		typedef mVector<int, R> 	iVector;

		public:
			Spline(unsigned bits_, F const &f_):
				loop(2),
				box(bits_), fd(f_)
			{}

			double operator()(fVector const &x)
			{
				return f(x);
			}

			double f(fVector const &x);
			double df(unsigned k, fVector const &x);
	};

	template <typename T>
	class Cache
	{
		std::queue<size_t> history;
		std::map<size_t, T> data;

		public:
			T &operator()(size_t i, std::function<T ()> const &f)
			{
				if (data.count(i) == 0)
				{
					if (history.size() > 65536)
					{
						size_t j = history.front();
						history.pop();
						data.erase(j);
					}

					history.push(i);
					data[i] = f();
				}

				return data[i];
			}
	};

	template <typename F>
	class Spline<F, 2>
	{
		enum { R = 2 };
		cVector<2> 	loop;
		cVector<2> 	b4;
		cVector<2> 	box;	
		F			fn;

		typedef mVector<double, R> 	fVector;
		typedef mVector<int, R> 	iVector;

		static double				A[256];
		mutable Cache<std::vector<double>>	cache;

		public:
			Spline(unsigned bits, F const &f_);
			double operator()(fVector const &x) const;
			double f(fVector const &x) const;
			double df(unsigned k, fVector const &x) const;
	};

	template <typename F>
	class Spline<F, 3>
	{
		enum { R = 3 };

		cVector<2> 	loop;
		cVector<R> 	b4;
		cVector<R> 	box;
		F		fn;
		
		typedef typename F::value_type value_type;
		typedef mVector<double, R> 	fVector;
		typedef mVector<int, R> 	iVector;

		static double			A[4096];
		mutable Cache<std::vector<value_type>>	cache;

		public:
			Spline(unsigned bits, F const &f_);
			value_type operator()(fVector const &x) const;
			value_type f(fVector const &x) const;
			value_type df(unsigned k, fVector const &x) const;
	};

#include "spline1.ih"
#include "spline2.ih"
#include "spline3.ih"

	} // namespace Interpol
} // namespace Misc

#include "cspline.ih"
#include "bspline.ih"
#include "tspline.ih"

// vim:sw=4:ts=4:tw=72
