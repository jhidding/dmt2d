#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>

#include "header.h"
#include "mregex.h"
#include "stauncharray.h"

using namespace std;
using namespace Conan;

Header::Header(istream &f)
{
	from_file(f);
}

inline string strip_newlines(string const &str)
{
	size_t pos = str.find_first_of('\n');
	return str.substr(0, pos);
}

void Header::to_file(ostream &f)
{
	ostringstream s;
	for (iterator i = begin(); i != end(); ++i)
	{
		s << i->first << "=" << strip_newlines(i->second) << endl;
	}

	string output(s.str());
	Conan::StaunchArray<char> raw(output.length());
	copy(output.begin(), output.end(), raw.begin());
	raw.to_file(f);
}

void Header::to_txt_file(ostream &f)
{
	f << "### Header ###" << endl;
	for (iterator i = begin(); i != end(); ++i)
	{
		f << "# " << i->first << "=" << strip_newlines(i->second) << endl;
	}
}

void Header::from_file(istream &f)
{
	clear();

	Conan::StaunchArray<char> raw(f);
	string input(raw.begin(), raw.end());

	string id_value_pattern = "([A-Za-z0-9_]*)=(.*)\n";
	Regex r(id_value_pattern, 10, REG_EXTENDED | REG_NEWLINE);

	while (input.length() > 0)
	{
		int pos = r.match(input);
		if (r.is_error())
		{
			cerr << "regex: " << r.error_msg() << endl;
			exit(0);
		}

		insert(value_type(r[1], r[2]));
		input = input.substr(pos);
		r.clear();
	}
}

