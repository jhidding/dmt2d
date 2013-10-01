/*
 * date() creates a date time stamp, for use in filenames
 */

#include "date.h"

#include <sstream>
#include <iomanip>
#include <ctime>

using namespace std;

string date()
{
	time_t tsec = time(NULL);
	struct tm *T = gmtime(&tsec);
	ostringstream s;
	s << setfill('0') 
		<< setw(2) << T->tm_year - 100
		<< setw(2) << T->tm_mon + 1 
		<< setw(2) << T->tm_mday; 
	return s.str();
}

