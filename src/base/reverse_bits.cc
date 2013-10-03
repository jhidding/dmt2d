#include "cvector.hh"
#include "misc.hh"

#include <algorithm>
#include <cstddef>

using namespace System;

ReverseBitorder::ReverseBitorder(unsigned b)
	: std::vector<size_t>(size_t(1) << b), bits(b)
{
	Range I(size());
	std::transform(I.begin(), I.end(), begin(), *this);
}

size_t ReverseBitorder::operator()(size_t i) const
{
	size_t r = 0, m = size_t(1) << bits;

	while (i)
	{
		m >>= 1;
		if (i & 1) r |= m;
		i >>= 1;
	}
	
	return r;
}

