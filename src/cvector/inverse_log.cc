#include "bit_tricks.h"

#include <algorithm>

using namespace BitTricks;


InverseLog::InverseLog(unsigned b)
	: std::vector<std::size_t>(std::size_t(1) << b), bits(b)
{
	std::generate(begin(), end(), Counter<std::size_t>(0));
	std::transform(begin(), end(), begin(), *this);
}

std::size_t InverseLog::operator()(std::size_t i) const
{
	std::size_t n = std::size_t(1) << bits;
	
	do {
		i >>= 1;
		n >>= 1;
	} while (i > 0);

	return n;
}

