#include "history.h"
#include "stauncharray.h"
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iomanip>

#include "convert.h"
#include "mregex.h"
#include "getopt.h"

using namespace std;
using namespace Conan;

inline string strip_newlines(string const &str)
{
	size_t pos = str.find_first_of('\n');
	return str.substr(0, pos);
}

History::History(istream &in)
{
	from_file(in);
	update();
}

void History::update()
{
	string command = GetOpt::instance().command_line();
	insert(pair<time_t, string>(time(NULL), command));
}

void History::to_file(ostream &out)
{
	ostringstream s;
	for (iterator i = begin(); i != end(); ++i)
	{
		s << i->first << ":" << strip_newlines(i->second) << endl;
	}

	string output(s.str());
	StaunchArray<char> raw(output.length());
	copy(output.begin(), output.end(), raw.begin());

	raw.to_file(out);
}

void History::to_txt_file(ostream &out)
{
	out << "### History ###" << endl;
	for (iterator i = begin(); i != end(); ++i)
	{
		std::string t(ctime(&(i->first)));
		t.erase(t.size() - 1);
		out << "# " << t << ": " << strip_newlines(i->second) << endl;
	}
}

void History::from_file(istream &in)
{
	StaunchArray<char> raw(in);
	string input(raw.begin(), raw.end());

	string id_value_pattern = "([0-9]*):(.*)\n";
	Regex r(id_value_pattern, 10, REG_EXTENDED | REG_NEWLINE);

	while (input.length() > 0)
	{
		int pos = r.match(input);
		if (r.is_error())
		{
			cerr << "regex: " << r.error_msg() << endl;
			exit(0);
		}

		time_t time = from_string<time_t>(r[1]);
		insert(value_type(time, r[2]));
		input = input.substr(pos);
		r.clear();
	}
}

