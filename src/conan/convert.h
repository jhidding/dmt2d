#pragma once

#include <string>
#include <sstream>

namespace Conan {

template <typename T>
std::string to_string(T const &A)
{
	std::ostringstream s;
	s << A;
	return s.str();
}

template <typename T>
T from_string(std::string const &B)
{
	T A;
	std::istringstream s(B);
	s >> A;
	return A;
}

}

