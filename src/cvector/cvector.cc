#include "cvector.h"
#include "cvector_type.h"

#include <iostream>

cVector::InfoMap cVector::info_map;

cVector::cVector(unsigned r, unsigned b_, std::size_t X_):
	info(info_map[cVector_type(r, b_)]), X(X_), bits(b_), rank(r)
{}


std::size_t cVector::size() const
{
	return info->size;
}

cVector &cVector::operator+=(cVector const &o)
{
	std::size_t Y = 0;
	for (unsigned i = 0; i < rank; ++i)
	{
		Y |= (X + ((o.X & o.info->mask[i]) << (info->bits[i] - o.info->bits[i]))) & info->mask[i];
	}

	X = Y;
	return *this;
}

cVector &cVector::operator-=(cVector const &o)
{
	std::size_t Y = 0;
	for (unsigned i = 0; i < rank; ++i)
	{
		Y |= (X - ((o.X & o.info->mask[i]) << (info->bits[i] - o.info->bits[i]))) & info->mask[i];
	}

	X = Y;
	return *this;
}

cVector cVector::reverse() const
{
	std::size_t n = 0;

	for (unsigned i = 0; i < rank; ++i)
	{
		n |= (info->rev[(*this)[i]] << info->bits[i]);
	}

	return cVector(rank, bits, n);
}

cVector cVector::reverse_morton() const
{
	std::size_t n = 0;

	for (unsigned i = 0; i < rank; ++i)
	{
		n |= (info->mor[(*this)[i]] << (rank - i - 1));
	}

	return cVector(bits, rank, n);
}

cVector cVector::morton() const
{
	std::size_t n = 0;

	for (unsigned i = 0; i < rank; ++i)
	{
		n |= (info->mor[(*this)[i]] << i);
	}

	return cVector(bits, rank, n);
}

cVector cVector::invlog() const
{
	std::size_t n = 0;

	for (unsigned i = 0; i < rank; ++i)
	{
		n |= (info->ilg[(*this)[i]] << info->bits[i]);
	}

	return cVector(rank, bits, n);
}

