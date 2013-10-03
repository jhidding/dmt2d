#include <iostream>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "../base/argv.hh"
#include "../base/array.hh"
#include "../base/boxconfig.hh"
#include "../base/date.hh"
#include "../base/header.hh"
#include "../base/history.hh"
#include "../base/misc.hh"
#include "../base/global.hh"

#include "../misc/hessian.hh"

using namespace System;

class LagrangianCatastropheData
{
	typedef mVector<double, 2> Point;
	Array<Point> umbilics;

	public:
		LagrangianCatastropheData(Array<Point> umbilics_):
			umbilics(umbilics_) {}

		void to_txt_file(std::ostream &out)
		{
			out << "# umbilics\n";
			for (auto p : umbilics) out << p << std::endl;
		}
};

Array<mVector<double, 2>> compute_double_zeros_2(ptr<BoxConfig<2>> box, Array<double> f, Array<double> g);

ptr<LagrangianCatastropheData> compute_lagrangian_catastrophes(ptr<BoxConfig<2>> box, Array<double> phi)
{
	DMT::Hessian<2> H(box, phi);
	size_t size = box->size();
	Array<double> Q1(size), Q2 = H[1];
	transform(H[0], H[2], Q1, [] (double a, double b) { return a - b; });

	std::vector<Array<double>> Eval(2, size);
	std::vector<Array<mVector<double,2>>> Evec(2, size);

	for (size_t x = 0; x < size; ++x)
	{
		double d = (H[0][x] + H[2][x]) / 2.0;
		double q = sqrt(pow(H[0][x] - H[2][x], 2) + 4 * pow(H[1][x], 2)) / 2.0;
		Eval[0][x] = d + q;
		Eval[1][x] = d - q;
		Evec[0][x][0] = Eval[0][x] - H[2][x]; 
		Evec[0][x][1] = H[1][x];
		Evec[1][x][0] = Eval[1][x] - H[2][x];
		Evec[1][x][1] = H[1][x];
	}

	auto umbilics = compute_double_zeros_2(box, Q1, Q2);

	return make_ptr<LagrangianCatastropheData>(umbilics);
}

void command_dmt(int argc_, char **argv_)
{
	Argv argv = read_arguments(argc_, argv_,
		Option({0, "h", "help", "false", "print this help."}),
		Option({Option::VALUED | Option::CHECK, "I", "id", date_string(),
			"Identifier for this run."}),
		Option({Option::VALUED | Option::CHECK, "i", "input", "-",
			"Input file. Defaults to <id>.init.dmt"}),
		Option({Option::VALUED | Option::CHECK, "o", "output", "-",
			"Output file. Defaults to <id>.catastrophes.init.dmt"}));

	std::string fn_input = (argv["input"] == "-" ? 
			timed_filename(argv["id"], "density", -1) :
			argv["input"]);

	std::string fn_output = (argv["output"] == "-" ?
			timed_filename(argv["id"], "catastrophes", -1) :
			argv["output"]);

	std::ifstream fi;
	std::cerr << "reading " << fn_input << " ...\n";
	fi.open(fn_input.c_str(), std::ios::in | std::ios::binary);

	System::Header 		H(fi);
	System::History 	I(fi);
	I << argv;

	unsigned	bits	= H.get<unsigned>("mbits");
	float		L 	= H.get<float>("size");

	auto box = make_ptr<BoxConfig<2>>(bits, L);

	std::ofstream fo;
	fo.open(fn_output.c_str(), std::ios::out | std::ios::binary);

	H.to_txt_file(fo);
	I.to_txt_file(fo);

	Array<double> potential = load_from_file<double>(fi, "potential");
	std::cerr << "box: " << box->N() << " data: " << potential.size() << std::endl;
	auto data = compute_lagrangian_catastrophes(box, potential);
	data->to_txt_file(fo);

	fi.close();
	fo.close();
}

System::Global<System::Command> COMMAND_DMT("lagrangian", command_dmt);

