#pragma once

#include <algorithm>
#include <map>
#include <queue>

#include "../conan/mvector.h"
#include "../conan/boxconfig.h"
#include "../cvector/cvector.h"

namespace Spline {

inline int integer(double a)
{
	if (a > 0)
		return static_cast<int>(a);
	else
		return static_cast<int>(a - 1);
}

inline double fraction(double a)
{
	return a - integer(a);
}

template <int rank>
inline Conan::mVector<int, rank> c2m(cVector c)
{
	Conan::mVector<int, rank> v;
	for (unsigned k = 0; k < rank; ++k) v[k] = c[k];
	return v;
}

template <int rank>
class Interpol
{
};

template <>
class Interpol<1>
{
	enum { rank = 1 };
	static double	A[16];
	unsigned	bits;
	double const	*fd;

	public:
		Interpol(double const *f_):
			bits(Conan::BoxConfig::bits()), fd(f_)
		{}

		double operator()(Conan::mVector<double, rank> const &x)
		{
			return f(x);
		}

		double f(Conan::mVector<double, rank> const &x)
		{
			Conan::mVector<double, rank> P;
			std::transform(x.begin(), x.end(), P.begin(), fraction);

			Conan::mVector<int, rank> start;
			std::transform(x.begin(), x.end(), start.begin(), integer);
			cVector S(rank, bits, Conan::BoxConfig::M(start));

			cVector dx(1, 1, 1);
			double b[4];
			for (cVector i(1, 1, 0); i < 2; ++i)
			{
				cVector X = S + i;
				b[i]   = fd[X];
				b[i+2] = (fd[X + dx] - fd[X - dx])/2;
			}

			double alpha[4];
			std::fill_n(alpha, 4, 0.0);

			for (cVector i(2, 2, 0); i < 16; ++i)
				alpha[i[1]] += b[i[0]] * A[i];

			double v = 0.0;
			for (unsigned i = 0; i < 4; ++i)
				v += alpha[i] * pow(P[0], i);

			return v;
		}

		double df(unsigned k, Conan::mVector<double, rank> const &x)
		{
			Conan::mVector<double, rank> P;
			std::transform(x.begin(), x.end(), P.begin(), fraction);

			Conan::mVector<int, rank> start;
			std::transform(x.begin(), x.end(), start.begin(), integer);
			cVector S(rank, bits, Conan::BoxConfig::M(start));

			cVector dx(1, 1, 1);
			double b[4];
			for (cVector i(1, 1, 0); i < 2; ++i)
			{
				cVector X = S + i;
				b[i]   = fd[X];
				b[i+2] = (fd[X + dx] - fd[X - dx])/2;
			}

			double alpha[4];
			std::fill_n(alpha, 4, 0.0);

			for (cVector i(2, 2, 0); i < 16; ++i)
				alpha[i[1]] += b[i[0]] * A[i];

			double v = 0.0;
			for (unsigned i = 1; i < 4; ++i)
				v += alpha[i] * i * pow(P[i], i - 1);

			return v;
		}
};

template <>
class Interpol<2>
{
	enum { rank = 2 };
	unsigned	bits;

	class Cache
	{
		static double			A[256];
		std::map<size_t, double *>	C;
		std::queue<size_t>		Q;
		double const			*f;
		size_t				cache_size;
		double				*cache_data;
		double				*cache_ptr;

		public:
			Cache(Cache const &o);
			Cache(double const *f_);
			double *operator[](cVector const &idx);
			void push(cVector const &S);
			double *pop();
			~Cache();
	} cache;

	public:
		Interpol(double const *f_);
		double operator()(Conan::mVector<double, rank> const &x);
		double f(Conan::mVector<double, rank> const &x);
		double df(unsigned k, Conan::mVector<double, rank> const &x);
};

template <>
class Interpol<3>
{
	enum { rank = 3 };
	unsigned		bits;
	
	class Cache
	{
		static double			A[4096];
		std::map<size_t, double *>	C;
		std::queue<size_t>		Q;
		double const			*f;
		size_t				cache_size;
		double				*cache_data;
		double				*cache_ptr;

		public:
			Cache(Cache const &o);
			Cache(double const *f_);
			double *operator[](cVector const &idx);
			void push(cVector const &S);
			double *pop();
			~Cache();
	} cache;

	public:
		Interpol(double const *f_);
		double operator()(Conan::mVector<double, rank> const &x);
		double f(Conan::mVector<double, rank> const &x);
		double df(unsigned k, Conan::mVector<double, rank> const &x);
};

} // namespace JH

