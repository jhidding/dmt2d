/*
 * DFT
 * very simple class wrapper for the 
 * wonderfull FFTW3 library
 */

#pragma once
#include <fftw3.h>
#include <complex>

class DFT
{
	std::complex<double>	*d_in, *d_out;
	fftw_plan		d_plan_fwd, d_plan_bwd;
	int			*d_dim;

	public:
		DFT(int rank, unsigned N);
		~DFT();
		void forward();
		void backward();
		std::complex<double> *in();
		std::complex<double> *out();
};

