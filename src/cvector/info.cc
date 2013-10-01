#include "cvector.h"

cVector::Info::Info(unsigned rank, unsigned b)
	: rev(b), mor(b, rank), ilg(b), size(std::size_t(1) << (b * rank))
{
	//std::cerr << "! running cvector_info !" << std::endl;
	bits = new unsigned[rank];
	mask = new std::size_t[rank];

	bits[0] = 0;
	mask[0] = (std::size_t(1) << b) - 1;
	for (unsigned i = 1; i < rank; ++i)
	{
		bits[i] = bits[i-1] + b;
		mask[i] = mask[i-1] << b;
	}
}

cVector::Info::~Info()
{
	delete[] bits;
	delete[] mask;
}

