#include "bit_tricks.h"

#include <algorithm>

using namespace BitTricks;

ReverseBitorder::ReverseBitorder(unsigned b)
	: std::vector<std::size_t>(std::size_t(1) << b), bits(b)
{
	std::generate(begin(), end(), Counter<std::size_t>(0));
	std::transform(begin(), end(), begin(), *this);
}

std::size_t ReverseBitorder::operator()(std::size_t i) const
{
	std::size_t r = 0, m = std::size_t(1) << bits;

	while (i)
	{
		m >>= 1;
		if (i & 1) r |= m;
		i >>= 1;
	}
	
	return r;
}

