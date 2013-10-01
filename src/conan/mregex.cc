#include "mregex.h"

Regex::Regex(std::string const &re, int mm_, int flags)
	: std::vector<std::string>(), max_match(mm_)
{
	erno = regcomp(&buffer, re.c_str(), flags);
	re_match = new regmatch_t[max_match];
}

int Regex::match(std::string const &s, int flags)
{
	erno = regexec(&buffer, s.c_str(), max_match, re_match, flags);

	if (!is_error())
	{
		for (int i = 0; i < max_match; ++i)
		{
			if (re_match[i].rm_so == -1) break;
			push_back(s.substr(
				re_match[i].rm_so,
				re_match[i].rm_eo - re_match[i].rm_so
			));
		}
		return re_match[0].rm_eo;
	}
	else return s.length();
}

bool Regex::is_error() { return erno != 0; }

std::string Regex::error_msg()
{
	regerror(erno, &buffer, erstr, 256);
	return std::string(erstr);
}

Regex::~Regex()
{
	delete[] re_match;
	regfree(&buffer);
}

