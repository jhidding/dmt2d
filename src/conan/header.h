/*
 * Reads one block from a file and interprets it as a header.
 * a header contains <name>=<value> pairs in ASCII separated
 * by semicolons ';'. The entire block is saved as a Fortran
 * array, encapsulated by a 32bit unsigned integer containing
 * the size of the block.
 */

#pragma once

#include <iosfwd>
#include <string>
#include <map>

namespace Conan
{

	class Header: public std::map<std::string, std::string>
	{
		public:
			Header() {}
			Header(std::istream &f);
			void to_file(std::ostream &f);
			void to_txt_file(std::ostream &f);
			void from_file(std::istream &f);
	};

} // namespace Conan
