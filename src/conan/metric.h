#pragma once

#include "mvector.h"
#include "modulus.h"
#include "boxconfig.h"

#include <vector>

namespace Conan
{

/*! 
 * \brief Calculates the position in an array
 * using quick bitshift operations
 */
class Metric
{
	unsigned N;
	unsigned rank;
	std::vector<unsigned> bits;

	public:
		Metric() {}

		Metric(unsigned rank_, unsigned b):
			N(1 << b), rank(rank_), bits(rank)
		{
			for (unsigned i = 0; i < rank; ++i)
				bits[i] = b * i;
		}

		template <int Rank>
		size_t operator()(mVector<int, Rank> const &x) const
		{
			size_t s = 0;
			for (unsigned i = 0; i < Rank; ++i)
				s |= modulus(x[i], N) << bits[i];
			return s;
		}

		template <int Rank>
		size_t operator()(mVector<double, Rank> const &x) const
		{
			size_t s = 0;
			for (unsigned i = 0; i < Rank; ++i)
				s |= modulus(static_cast<int>(floor(x[i])), N) << bits[i];
			return s;
		}
};

} // namespace Conan

