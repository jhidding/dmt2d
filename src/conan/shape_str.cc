#include "shape_str.h"

#include <sstream>

std::string Conan::shape_str(size_t N, size_t rank)
{
	std::ostringstream s;
	s << "(";
	for (size_t k = 0; k < rank - 1; ++k)
		s << N << ", ";
	s << N << ")";
	return s.str();
}

