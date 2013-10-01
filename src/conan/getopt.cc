#include "getopt.h"

#include <sstream>
#include <iterator>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace Conan;

GetOpt *GetOpt::s_singleton = NULL;

GetOpt::GetOpt(int argc_, char **argv_)
	: argc(argc_), argv(argv_)
{
	copy(argv, argv + argc, back_inserter(arguments));
}

GetOpt &GetOpt::initialize(int argc_, char **argv_)
{
	if (s_singleton == NULL)
	{
		s_singleton = new GetOpt(argc_, argv_);
	}

	return *s_singleton;
}

string GetOpt::command_line()
{
	ostringstream S;
	copy(arguments.begin(), arguments.end(), ostream_iterator<string>(S, " "));
	return S.str();
}

string GetOpt::name()
{
	return arguments[0];
}

size_t GetOpt::size()
{
	return arguments.size();
}

string GetOpt::operator[](size_t idx)
{
	return arguments[idx];
}

string GetOpt::get(std::string const &s)
{
	vector<string>::iterator loc = find(arguments.begin(), arguments.end(), s);
	if (loc >= arguments.end() - 1) return string("");
	else return *(++loc);
}

bool GetOpt::has(std::string const &s)
{
	return find(arguments.begin(), arguments.end(), s) != arguments.end();
}

