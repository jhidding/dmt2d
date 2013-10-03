#include "cvector.hh"

#include <algorithm>
#include <cstddef>

using namespace System;

InverseLog::InverseLog(unsigned b)
	: std::vector<size_t>(size_t(1) << b), bits(b)
{
	Range I(size());
	std::transform(I.begin(), I.end(), begin(), *this);
}

size_t InverseLog::operator()(size_t i) const
{
	size_t n = size_t(1) << bits;
	
	do {
		i >>= 1;
		n >>= 1;
	} while (i > 0);

	return n;
}

