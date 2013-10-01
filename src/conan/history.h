#pragma once

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <map>

namespace Conan
{
	class History: public std::multimap<time_t, std::string>
	{
		public:
			History() {}
			History(std::istream &in);
			void update();
			void from_file(std::istream &in);
			void to_file(std::ostream &out);
			void to_txt_file(std::ostream &out);
	};
}

