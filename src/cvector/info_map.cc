#include "cvector.h"
#include "cvector_type.h"

cVector::Info const *cVector::InfoMap::operator[](unsigned idx)
{
	#pragma omp critical
	{
		if (count(idx) == 0)
		{
			cVector_type_union U;
			U.value = idx;

			insert(std::pair<unsigned, Info *>(idx, new Info(U.p.rank, U.p.bits)));
		}
	}

	return std::map<unsigned, Info *>::operator[](idx);
}

cVector::InfoMap::~InfoMap()
{
	for (iterator i = begin(); i != end(); ++i)
	{
		delete i->second;
	}
}

