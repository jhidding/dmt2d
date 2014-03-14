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
#include "../misc/interpol.hh"

using namespace System;

Array<mVector<double, 2>> compute_double_zeros_2(ptr<BoxConfig<2>> box, Array<double> f, Array<double> g);

Array<double> LinSpace(double a, double b, unsigned n)
{
	Array<double> r(n);
	for (unsigned i = 0; i < n; ++i)
		r[i] = a + ((b - a) * i) / (n - 1);
	return r;
}

template <typename A>
typename A::value_type mean(A const &a)
{
	return std::accumulate(a.begin(), a.end(), typename A::value_type(0)) / a.size();
}

void compute_pdf(ptr<BoxConfig<2>> box, Array<double> rho, 
	Array<double> phi, double aa, unsigned n, std::ostream &fo)
{
	std::cerr << "hessian ... ";
	DMT::Hessian<2> H(box, phi);

	std::cerr << "sigma: " << mean(map(rho, [] (double a) { return a*a; })) << std::endl;

	size_t size = box->size();
	Array<double> Q1(size), Q2 = H[1];
	transform(H[0], H[2], Q1, [] (double a, double b) { return a - b; });

	cVector<2> b(box->bits());
	for (size_t x = 0; x < size; ++x)
	{
		double p = (H[0][x] + H[2][x]) / 2.0;
		double q = sqrt(pow(H[0][x] - H[2][x], 2) + 4 * pow(H[1][x], 2)) / 2.0;
		fo << b.dvec(x) << " " << rho[x] << " " << p << " " << q << std::endl;
	}
	
	fo << "\n\n\n";

	Array<double> A_level(size), d_level(size);
	for (unsigned i = 0; i < n; ++i)
	{
		double A = -aa + (2*aa*i)/(n-1);
		for (size_t x = 0; x < size; ++x)
			A_level[x] = (H[0][x] + H[2][x]) / 2.0 - A;

		for (unsigned j = 1; j < n; ++j)
		{
			double d = (aa * j)/(n-1);
			for (size_t x = 0; x < size; ++x)
				d_level[x] = sqrt(pow(H[0][x] - H[2][x], 2) + 4 * pow(H[1][x], 2)) / 2.0 - d;

			auto pts = compute_double_zeros_2(box, A_level, d_level);
			fo << A-d << " " << A+d << " " << pts.size() / (box->L() * box->L()) << " " 
				<< A << " " << d << std::endl;
		}
		fo << std::endl;
	}
	//auto umbilics = compute_double_zeros_2(box, Q1, Q2);
}

void command_pdf(int argc_, char **argv_)
{
	Argv argv = read_arguments(argc_, argv_,
		Option({0, "h", "help", "false", "print this help."}),
		Option({Option::VALUED | Option::CHECK, "I", "id", date_string(),
			"Identifier for this run."}),
		Option({Option::VALUED | Option::CHECK, "i", "input", "-",
			"Input file. Defaults to <id>.init.dmt"}),
		Option({Option::VALUED | Option::CHECK, "a", "bound", "1.0",
			"maximum deviation value for eigenvalues."}),
		Option({Option::VALUED | Option::CHECK, "n", "bins", "51",
			"number of bins for eigenvalues."}),
		Option({Option::VALUED | Option::CHECK, "o", "output", "-",
			"Output file. Defaults to <id>.pdf.init.dmt"}));

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
	double 		a	= argv.get<double>("bound");
	unsigned 	n 	= argv.get<unsigned>("bins");

	auto box = make_ptr<BoxConfig<2>>(bits, L);

	std::ofstream fo;
	fo.open(fn_output.c_str(), std::ios::out | std::ios::binary);

	H.to_txt_file(fo);
	I.to_txt_file(fo);

	Array<double> 
		density = load_from_file<double>(fi, "density"),
		potential = load_from_file<double>(fi, "potential");

	std::cerr << "box: " << box->N() << " data: " << potential.size() << std::endl;
	compute_pdf(box, density, potential, a, n, fo);
	std::cerr << "writing to file ... ";

	fi.close();
	fo.close();
}

System::Global<System::Command> COMMAND_PDF("pdf", command_pdf);

