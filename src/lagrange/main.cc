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
Array<Array<mVector<double, 2>>> compute_level_set_2(ptr<BoxConfig<2>> box, double a, Array<double> f);
void print_level_set_2(ptr<BoxConfig<2>>, std::ostream &, Array<Array<mVector<double, 2>>>);
Array<Array<mVector<double, 2>>> compute_a3_lines_2(ptr<BoxConfig<2>>, Array<double>, DMT::Hessian<2>, Array<mVector<double,2>>);
Array<mVector<double, 2>> compute_line_extrema_2(ptr<BoxConfig<2>>, Array<double>, Array<Array<mVector<double,2>>> a3);
Array<Array<mVector<double, 2>>> compute_pogosyan_lines_2(ptr<BoxConfig<2>>, Array<double>, Array<mVector<double,2>>, Array<mVector<double,2>>, DMT::Hessian<2>, bool);

template <typename F>
void print_a3_line_2(ptr<BoxConfig<2>> box, std::ostream &out, Array<mVector<double,2>> C, F f)
{
	bool pen_on_paper = true;
	for (unsigned a = 0; a < C.size(); ++a)
	{
		unsigned i = System::modulus(int(a), int(C.size())),
			 j = System::modulus(int(a) - 1, int(C.size()));

		if ((C[i] - C[j]).norm() > box->L() / 2)
			out << "\n\n";

		if (f(C[i]) > 0)
		{
			pen_on_paper = true;
			out << C[i] << " " << f(C[i]) << "\n";
		}
		else
		{
			if (pen_on_paper)
				out << "\n\n";
			pen_on_paper = false;
		}
	}
}

template <typename F>
void print_entire_a3_line_2(ptr<BoxConfig<2>> box, std::ostream &out, Array<mVector<double,2>> C, F f)
{
	for (unsigned a = 0; a < C.size(); ++a)
	{
		unsigned i = System::modulus(int(a), int(C.size())),
			 j = System::modulus(int(a) - 1, int(C.size()));

		if ((C[i] - C[j]).norm() > box->L() / 2)
			out << "\n\n";

		out << C[i] << " " << f(C[i]) << "\n";
	}
}


class LagrangianCatastropheData
{
	ptr<BoxConfig<2>> box;
	typedef mVector<double, 2> Point;
	typedef Array<Point> Line;

	typedef Array<Point> Point_set;
	typedef Array<Line> Line_set;

	Array<double> rho, ev_a, ev_b;
	Array<Point> umbilics;
	Line_set a3a, a3b, pogo_a, pogo_b;
	Point_set Ea, Eb;

	public:
		LagrangianCatastropheData(ptr<BoxConfig<2>> box_, Array<double> rho_,
				Array<double> ev_a_, Array<double> ev_b_,
				Line_set p1, Line_set p2, Line_set Q1_, Line_set Q2_,
				Point_set ea_, Point_set eb_, Point_set umbilics_):
			box(box_), rho(rho_), ev_a(ev_a_), ev_b(ev_b_), umbilics(umbilics_), 
			a3a(Q1_), a3b(Q2_), pogo_a(p1), pogo_b(p2),
       			Ea(ea_), Eb(eb_)
		{}

		void to_txt_file(std::ostream &out, std::string const &id)
		{
			Misc::Interpol::Linear<Array<double>,2> Eva(box, ev_a), Evb(box, ev_b), Rho(box, rho);

			cVector<2> b(box->bits());
			for (size_t x = 0; x < box->size(); ++x)
			{
				auto X = b.dvec(x);
				out << X << rho[x] << " " << ev_a[x] << " " << ev_b[x] << std::endl;
			}

			std::ofstream fo1(timed_filename(id, "points", -1));
			fo1 << "# umbilics\n";
			for (auto p : umbilics) if (Eva(p) > 0) fo1 << p << " " << Eva(p) << std::endl;

			fo1 << "\n\n# a3+/- and a4 for alpha\n";
			for (auto p : Ea) if (Eva(p) > 0) fo1 << p << std::endl;

			fo1 << "\n\n# a3+/- and a4 for beta\n";
			for (auto p : Eb) if (Evb(p) > 0) fo1 << p << std::endl;
			fo1.close();

			std::ofstream fo2(timed_filename(id, "a3.alpha", -1));
			for (auto L : a3a) 
			{
				print_a3_line_2(box, fo2, L, Eva);
				fo2 << "\n\n";
			} fo2.close();

			/*
			std::ofstream fo4(timed_filename(id, "pogo.a", -1));
			std::ofstream fo5(timed_filename(id, "pogo.b", -1));
			for (auto L : pogo_a)
			{
				print_a3_line_2(box, fo4, L, Rho);
				fo4 << "\n\n";
			}
			for (auto L : pogo_b)
			{
				print_a3_line_2(box, fo5, L, Rho);
				fo5 << "\n\n";
			}
			fo4.close(); fo5.close();
			*/

			std::ofstream fo3(timed_filename(id, "a3.beta", -1));
			for (auto L : a3b) 
			{
				print_a3_line_2(box, fo3, L, Evb);
				fo3 << "\n\n";
			} fo3.close();

		}
};

ptr<LagrangianCatastropheData> compute_lagrangian_catastrophes(ptr<BoxConfig<2>> box, Array<double> rho, Array<double> phi)
{
	std::cerr << "hessian ... ";
	DMT::Hessian<2> H(box, phi);
	DMT::Hessian<2> H_rho(box, rho);
	std::cerr << "q1, q2 ... ";
	size_t size = box->size();
	Array<double> Q1(size), Q2 = H[1];
	transform(H[0], H[2], Q1, [] (double a, double b) { return a - b; });

	std::cerr << "eigenvalues ... ";
	std::vector<Array<double>> Eval; 
	Eval.push_back(Array<double>(size)); Eval.push_back(Array<double>(size));
	std::vector<Array<double>> Eval_rho;
	Eval_rho.push_back(Array<double>(size)); Eval_rho.push_back(Array<double>(size));
	std::vector<Array<mVector<double,2>>> Evec;
	Evec.push_back(Array<mVector<double,2>>(size)); Evec.push_back(Array<mVector<double,2>>(size));

	for (size_t x = 0; x < size; ++x)
	{
		// phi
		double d = (H[0][x] + H[2][x]) / 2.0;
		double q = sqrt(pow(H[0][x] - H[2][x], 2) + 4 * pow(H[1][x], 2)) / 2.0;
		Eval[0][x] = d + q;
		Eval[1][x] = d - q;

		// rho
		d = (H_rho[0][x] + H_rho[2][x]) / 2.0;
		q = sqrt(pow(H_rho[0][x] - H_rho[2][x], 2) + 4 * pow(H_rho[1][x], 2)) / 2.0;
		Eval_rho[0][x] = d + q;
		Eval_rho[1][x] = d - q;
		Evec[0][x][0] = Eval_rho[0][x] - H_rho[2][x]; 
		Evec[0][x][1] = H_rho[1][x];
		Evec[1][x][0] = Eval_rho[1][x] - H_rho[2][x];
		Evec[1][x][1] = H_rho[1][x];
	}

	auto umbilics = compute_double_zeros_2(box, Q1, Q2);
	//auto Q1_zero = compute_level_set_2(box, 0.0, Q1);
	//auto Q2_zero = compute_level_set_2(box, 0.0, Q2);
	auto A3_alpha = compute_a3_lines_2(box, Eval[0], H, umbilics);
	auto A3_beta  = compute_a3_lines_2(box, Eval[1], H, umbilics);
	auto extr_alpha = compute_line_extrema_2(box, Eval[0], A3_alpha);
	auto extr_beta = compute_line_extrema_2(box, Eval[1], A3_beta);

	auto pogo_lines_a = compute_pogosyan_lines_2(box, rho, Evec[0], Evec[1], H_rho, true);
	auto pogo_lines_b = compute_pogosyan_lines_2(box, rho, Evec[0], Evec[1], H_rho, false);

	std::cerr << "[returning]\n";
	return make_ptr<LagrangianCatastropheData>(box, rho, Eval[0], Eval[1], pogo_lines_a, pogo_lines_b, A3_alpha, A3_beta, extr_alpha, extr_beta, umbilics);
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

	Array<double> 
		density = load_from_file<double>(fi, "density"),
		potential = load_from_file<double>(fi, "potential");

	std::cerr << "box: " << box->N() << " data: " << potential.size() << std::endl;
	auto data = compute_lagrangian_catastrophes(box, density, potential);
	std::cerr << "writing to file ... ";
	data->to_txt_file(fo, argv["id"]);

	fi.close();
	fo.close();
}

System::Global<System::Command> COMMAND_DMT("lagrangian", command_dmt);

