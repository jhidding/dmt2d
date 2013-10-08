#ifdef UNITTEST
#include "unittest.hh"
#include "cvector.hh"
#include <iostream>

using namespace System;

Test::Unit CVector_test("0000 - cvector",
	"CVector collects a set of functions to compute "
	"indices into power-2 data cubes.",
	[] ()
{
	cVector<2> b(3);

	for (size_t i = 0; i < b.size(); ++i)
	{
		size_t j = b.sub(i, b.sq_i[3]);
		size_t k = b.add(i, b.dx_i[0]);

		std::cerr << b.dvec(i) << " | " << b.dvec(j) << " | " << b.dvec(k) << std::endl;
	}

	return true;
});

#endif

