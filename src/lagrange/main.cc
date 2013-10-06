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
	ptr<BoxConfig<2>> box;
	typedef mVector<double, 2> Point;
	typedef Array<Point> Contour;
	typedef Array<Contour> Level_set;

	Array<Point> umbilics;
	Level_set Q1z, Q2z;

	public:
		LagrangianCatastropheData(ptr<BoxConfig<2>> box_, Level_set Q1_, Level_set Q2_, Array<Point> umbilics_):
			box(box_), umbilics(umbilics_), Q1z(Q1_), Q2z(Q2_) {}

		void to_txt_file(std::string const &id)
		{
			std::ofstream fo1(timed_filename(id, "umbilics", -1));
			out << "# umbilics\n";
			for (auto p : umbilics) fo1 << p << std::endl;
			fo1.close();

			std::ofstream fo2(timed_filename(id, "q1", -1));
			print_level_set_2(box, fo2, Q1z); fo2.close();

			std::ofstream fo3(timed_filename(id, "q2", -1));
			print_level_set_2(box, fo3, Q2z); fo3.close();
		}
};

Array<mVector<double, 2>> compute_double_zeros_2(ptr<BoxConfig<2>> box, Array<double> f, Array<double> g);
Array<Array<mVector<double, 2>>> compute_level_set_2(ptr<BoxConfig<2>> box, double a, Array<double> f);

ptr<LagrangianCatastropheData> compute_lagrangian_catastrophes(ptr<BoxConfig<2>> box, Array<double> phi)
{
	std::cerr << "hessian ... ";
	DMT::Hessian<2> H(box, phi);
	std::cerr << "q1, q2 ... ";
	size_t size = box->size();
	Array<double> Q1(size), Q2 = H[1];
	transform(H[0], H[2], Q1, [] (double a, double b) { return a - b; });

	std::cerr << "eigenvalues ... ";
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

	std::cerr << "umbilics ... ";
	auto umbilics = compute_double_zeros_2(box, Q1, Q2);
	auto Q1_zero = compute_level_set_2(box, 0.0, Q1);
	auto Q2_zero = compute_level_set_2(box, 0.0, Q2);

	std::cerr << "[returning]\n";
	return make_ptr<LagrangianCatastropheData>(box, Q1_zero, Q2_zero, umbilics);
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
	std::cerr << "writing to file ... ";
	data->to_txt_file(argv["id"]);

	fi.close();
	fo.close();
}

System::Global<System::Command> COMMAND_DMT("lagrangian", command_dmt);

