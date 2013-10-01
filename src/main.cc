#include <iostream>
#include <cstdlib>

#include "conan/boxconfig.h"
#include "conan/getopt.h"
#include "conan/history.h"
#include "conan/default_filename.h"
#include "dmt.h"

using namespace Conan;
using namespace std;

void usage()
{
	cout << "usage: dmt -input <filename> -output <filename>\n";
}

int main(int argc, char **argv)
{
	GetOpt &options = GetOpt::initialize(argc, argv);

	if (options.has("-help"))
	{
		usage();
		exit(0);
	}

	string fn_input = options.get_value("-input", timed_filename("delta0"));
	string fn_output = options.get_value("-output", default_filename("dmt"));

	ifstream fi;
	fi.open(fn_input.c_str(), ios::in | ios::binary);

	Header 		H(fi);
	History 	I(fi);

	unsigned	rank 	= from_string<unsigned>(H["rank"]);
	unsigned	N	= from_string<unsigned>(H["N"]);
	float		L 	= from_string<float>(H["L"]);

	BoxConfig::initialize(rank, N, L);

	ofstream fo;
	fo.open(fn_output.c_str(), ios::out | ios::binary);

	H.to_txt_file(fo);
	I.to_txt_file(fo);

	try {
		switch (rank)
		{
			//case 1:	nerve<1>(fi, fo); break;
			case 2:	dmt<2>(fi, fo); break;
			//case 3:	dmt<3>(fi, fo); break;
		}
	}
	catch (char const *msg)
	{
		cerr << msg << endl;
		exit(1);
	}

	fi.close();
	fo.close();

	return 0;
}

