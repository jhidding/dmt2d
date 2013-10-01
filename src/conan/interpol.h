#pragma once
#include "../cvector/cvector.h"
#include "boxconfig.h"
#include "mvector.h"

namespace Conan {
inline double ci_fraction(double a)
{
	return a - floor(a);
}

inline int ci_integer(double a)
{
	return static_cast<int>(floor(a));
}

template <typename T, int rank>
class Interpol
{
	T const		*f;
	Conan::Metric		M;

	public:
		Interpol() {}
		Interpol(T const *f_):
			f(f_)
		{}

		T operator()(mVector<double, rank> const &x) const
		{
			mVector<double, rank> A[2];
			std::transform(x.begin(), x.end(), A[1].begin(), ci_fraction);
			std::fill(A[0].begin(), A[0].end(), 1.0);
			A[0] -= A[1];

			mVector<int, rank> start;
			std::transform(x.begin(), x.end(), start.begin(), ci_integer);
			cVector S(rank, Conan::BoxConfig::bits(), Conan::BoxConfig::M(start));

			T v = 0;
			for (cVector i(rank, 1, 0); i < (1 << rank); ++i)
			{
				float z = 1;
				for (unsigned k = 0; k < rank; ++k) z *= A[i[k]][k];
				v += f[S + i] * z;
			}

			return v;
		}
};
}
