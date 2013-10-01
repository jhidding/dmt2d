#pragma once

#include <string>
#include <vector>
#include <sstream>

#include "convert.h"

namespace Conan {

class GetOpt
{
	static GetOpt *s_singleton;	
	std::vector<std::string> arguments;

	public:
		int argc;
		char **argv;

		static GetOpt &initialize(int argc_, char * argv_[]);
		static GetOpt &instance() { return *s_singleton; }

		std::string command_line();
		std::string name();

		size_t size();
		std::string operator[](size_t idx);
		std::string get(std::string const &s);

		template <typename T> 
		T get_value(std::string const &s, T def);

		bool has(std::string const &s);

	private:
		GetOpt(int argc_, char **argv_);
};


template <typename T> 
T GetOpt::get_value(std::string const &s, T def)
{ 
	return (has(s) ? from_string<T>(get(s)) : def);
}

}


