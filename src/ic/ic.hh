#pragma once
#include "../base/system.hh"
#include <memory>

namespace Conan
{
	//	unsigned    mbits = C.get<unsigned>("mbits");
	//	unsigned      dim = C.get<unsigned>("dim");
	//	double	        L = C.get<double>("L");
	//	double 	    slope = C.get<double>("power-slope");
	//	bool       smooth = C.get<bool>("smooth");
	//	double      sigma = C.get<double>("scale");
	//	unsigned     seed = C.get<unsigned>("seed");
	extern System::Array<double> generate_random_field(System::Header const &C);
	extern void compute_potential(System::Header const &C, System::Array<double>);
	extern void compute_displacement(System::Header const &C, System::Array<double>,
		std::ostream &fo);
}

