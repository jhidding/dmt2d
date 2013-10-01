#include "default_filename.h"
#include "getopt.h"
#include "date.h"

#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace Conan;

std::string default_filename(std::string const &q)
{
	GetOpt options = GetOpt::instance();
	string id = options.get_value("-id", date());
	string fn = id + "." + q + ".conan";
	return fn;
}

std::string adhesion_filename(float b)
{
	GetOpt options = GetOpt::instance();
	string id = options.get_value("-id", date());
	ostringstream s;
	s 	<< id + ".wg." << setfill('0') << setw(4)
		<< static_cast<int>(round(b * 1000)) 
		<< ".conan";
	return s.str();
}

std::string timed_filename(std::string const &q, float b)
{
	GetOpt options = GetOpt::instance();
	string id = options.get_value("-id", date());

	ostringstream s;

	if (b < 0.0)
	{
		s 	<< id << "." << q << "." << "init" 
			<< ".conan";
	}
	else
	{
		s 	<< id << "." << q << "." << setfill('0') << setw(4)
			<< static_cast<int>(round(b * 1000)) 
			<< ".conan";
	}
	return s.str();
}

