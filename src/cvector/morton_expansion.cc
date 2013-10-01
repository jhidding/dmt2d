#include "bit_tricks.h"

#include <algorithm>

using namespace BitTricks;

MortonExpansion::MortonExpansion(unsigned r, unsigned b):
	std::vector<std::size_t>(std::size_t(1) << r), rank(r), bits(b)
{
	std::generate(begin(), end(), Counter<std::size_t>(0));
	std::transform(begin(), end(), begin(), *this);
}

std::size_t MortonExpansion::operator()(std::size_t x)
{
	std::size_t m = 0;
	for (unsigned i = 0; i < rank; ++i)
	{
		m |= (x & (std::size_t(1) << (i * bits)));
		x <<= (bits - 1);
	}
	return m;
}

