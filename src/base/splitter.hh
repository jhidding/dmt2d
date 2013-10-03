#pragma once
#include <string>
#include <iterator>
#include <iostream>

namespace Misc
{
	class Splitter
	{
		std::string data;
		char delim;

		public:
			class const_iterator;
			Splitter(std::string const &s, char delim_);
			const_iterator begin() const;
			const_iterator end() const;
	};

	class Splitter::const_iterator: 
		public std::iterator<std::forward_iterator_tag, std::string>
	{
		Splitter const &S;
		size_t p1, p2;

		public:
			const_iterator(Splitter const &S_, size_t p = std::string::npos);
			const_iterator &operator++();
			std::string operator*() const;
			bool operator==(const_iterator const &o) const;
			bool operator!=(const_iterator const &o) const;
	};
}

// vim:ts=4:sw=4:tw=80:fdm=marker
