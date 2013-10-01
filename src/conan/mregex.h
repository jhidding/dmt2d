#pragma once

#include <regex.h>
#include <sys/types.h>
#include <vector>
#include <string>

class Regex: public std::vector<std::string>
{
	regex_t buffer;
	regmatch_t *re_match;
	int erno, max_match;
	char erstr[256];

	public:
		Regex(std::string const &re, int mm_=10, int flags=REG_EXTENDED);
		int match(std::string const &s, int flags=0);
		bool is_error();
		std::string error_msg();
		~Regex();
};

