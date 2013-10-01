#pragma once
#include <iostream>
#include "operators.h"
#include "level_set.h"
#include "conan/default_filename.h"
#include "conan/gridspace.h"
#include "conan/interpol.h"
#include "cvector/cvector.h"
#include "conan/getopt.h"
#include "gradient.h"
#include <cmath>

void output_cube_for_gnuplot(std::ostream &out, double *f)
{
	size_t const 	size 	= BoxConfig::size();
	int const 	N 	= BoxConfig::N();

	for (size_t i = 0; i < size; ++i)
	{
		if (i % N == 0)
			out << "\n";

		out << f[i] << " ";
	}
}

template <int R>
void dmt(std::istream &in, std::ostream &out)
{
	Conan::GetOpt &options = Conan::GetOpt::instance();
	size_t const 	size 	= BoxConfig::size();
	int const 	N 	= BoxConfig::N();
	double const 	L 	= BoxConfig::L();
//	int const 	bits	= BoxConfig::bits();

	double time = options.get_value("-b", 1.0);
	double corr = options.get_value("-c", 1.0);

	DFT dft(R, N);				// fftw wrapper
	constexpr int SMS = (R * R + R) / 2;	// symmetric matrix size
	Conan::KSpace<R> K(N, L);		// K-space

	Conan::Cube<double> delta_0(in), Phi;
	Conan::Cube<double> H[SMS];

	complex64 *d_hat = new complex64[size];
	std::copy(delta_0.begin(), delta_0.end(), dft.in());
	dft.forward();
	std::copy(dft.out(), dft.out() + size, d_hat);
	calc_hessian<R>(dft, d_hat, H);
	calc_potential<R>(dft, d_hat, Phi);
	//calc_displacement<R>(dft, d_hat, psi);
	Gradient<R> psi(Phi);
	delete[] d_hat;

	Conan::Cube<double> Q1;
	std::transform(H[0].begin(), H[0].end(), H[2].begin(), 
			Q1.begin(), std::minus<double>());

	// get eigenvalues and eigenvectors
	Conan::Cube<double> Eval[2];
	Conan::Cube<double> Evec[2][2];
	for (size_t x = 0; x < size; ++x)
	{
		double d = (H[0][x] + H[2][x]) / 2.0;
		double q = sqrt(pow(H[0][x] - H[2][x], 2) + 4 * pow(H[1][x], 2)) / 2.0;
		Eval[0][x] = d + q;
		Eval[1][x] = d - q;
		Evec[0][0][x] = Eval[0][x] - H[2][x]; 
		Evec[0][1][x] = H[1][x];
		Evec[1][0][x] = Eval[1][x] - H[2][x];
		Evec[1][1][x] = H[1][x];
	}

	double *ev_minus_1 = new double[size];
	std::transform(Eval[0].begin(), Eval[0].end(), ev_minus_1,
		[time] (double a) { return a - 1./time; });
	Level_set_2 caustics0(ev_minus_1);

	std::transform(Eval[1].begin(), Eval[1].end(), ev_minus_1,
		[time] (double a) { return a - 1./time; });
	Level_set_2 caustics1(ev_minus_1);
	delete[] ev_minus_1;
	
	Conan::Interpol<mVector<double, 2>, R> vel(psi.data());
	std::function<mVector<double, 2> (mVector<double, 2> const &)>
		zeldovich = [time, corr, L, N, &vel] (mVector<double, 2> const &p)
	{
		return (p * BoxConfig::scale() + vel(p) * time * corr) % BoxConfig::L();
	};

	Level_set_2 zero_q1(Q1.data()), zero_q2(H[1].data());
	Double_level_set_2 zero_q12(Q1.data(), H[1].data());
	A3_lines_2 
		a3_a(Eval[0].data(), H[0].data(), H[1].data(), 
			H[2].data(), zero_q12, true, true),
		a3_b(Eval[1].data(), H[0].data(), H[1].data(), 
			H[2].data(), zero_q12, true, true);

	std::vector<mVector<double, 2>> a4 = a3_a.A4_points(Eval[0].data());

	Conan::GridSpace<double, R> G(N, L), NG(N, N);
	for (size_t x = 0; x < size; ++x)
	{
		out << G[x] << " " << zeldovich(NG[x]) << " "
			<< Eval[0][x] << " " << Eval[1][x] << "\n";
	}	
	out << "\n\n\n";
	zero_q12.output_for_gnuplot(Eval[0].data(), out);
	out << "\n\n\n";

	std::for_each(a4.begin(), a4.end(),
		[&] (mVector<double, 2> const &x)
	{
		size_t i = BoxConfig::M(x);
		if (Eval[0][i] > 0)
			out << x << " " << Eval[0][i] << "\n";
	});
	out << "\n\n\n";
	zero_q12.map(zeldovich).output_for_gnuplot(out);
	out << "\n\n\n";
	output_cube_for_gnuplot(out, Eval[0].data());

	// contour output must go to separate files, due to
	// gnuplot limitations
	std::ofstream 
		of_cntr1(default_filename("cntr1.dmt")),
		of_cntr2(default_filename("cntr2.dmt")),
		of_a3aL(default_filename("a3a.L.dmt")),
		of_a3bL(default_filename("a3b.L.dmt")),
		of_a3aE(default_filename("a3a.E.dmt")),
		of_a3bE(default_filename("a3b.E.dmt")),
		of_caustic0(default_filename("caustic0.dmt")),
		of_caustic1(default_filename("caustic1.dmt"));

	caustics0.output_for_gnuplot(of_cntr1);
	caustics1.output_for_gnuplot(of_cntr2);
	a3_a.output_for_gnuplot(Eval[0].data(), of_a3aL);
	a3_b.output_for_gnuplot(Eval[1].data(), of_a3bL);
	a3_a.output_mapped(zeldovich, Eval[0].data(), of_a3aE);
	a3_b.output_mapped(zeldovich, Eval[1].data(), of_a3bE);

	std::for_each(caustics0.begin(), caustics0.end(),
		[&] (Contour_2 const &C)
	{
		C.map(zeldovich).output_for_gnuplot(of_caustic0);
		of_caustic0 << "\n\n";
	});
	std::for_each(caustics1.begin(), caustics1.end(),
		[&] (Contour_2 const &C)
	{
		C.map(zeldovich).output_for_gnuplot(of_caustic1);
		of_caustic1 << "\n\n";
	});

	std::ofstream of_evec(default_filename("evec.dmt"));
	for (cVector x(2, BoxConfig::bits(), 0); x < size; ++x)
	{
		if (x[0] % 4 != 0 or x[1] % 4 != 0) continue;

		of_evec << x[0] << " " << x[1] << " " << Eval[0][x] - H[2][x] << " "
			<< H[1][x] << " " << Eval[0][x] - H[0][x] << "\n";
	}
	of_evec.close();

	//of_cntr1.close(); of_cntr2.close(); of_a3a.close(); of_a3b.close();
}

