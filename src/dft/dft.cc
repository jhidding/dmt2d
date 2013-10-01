#include "dft.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>

DFT::DFT(int rank, unsigned N)
{
	size_t size = 1;
	for (int i = 0; i < rank; ++i) size *= N;
	size_t byte_size = size * sizeof(std::complex<double>);

	d_dim = new int[rank];
	std::fill_n(d_dim, rank, N);

	d_in = reinterpret_cast<std::complex<double> *>(fftw_malloc(byte_size));
	d_out = reinterpret_cast<std::complex<double> *>(fftw_malloc(byte_size));

	d_plan_fwd = fftw_plan_dft(rank, d_dim,
		reinterpret_cast<fftw_complex *>(d_in), 
		reinterpret_cast<fftw_complex *>(d_out), FFTW_FORWARD, FFTW_ESTIMATE);

	d_plan_bwd = fftw_plan_dft(rank, d_dim,
		reinterpret_cast<fftw_complex *>(d_in),
		reinterpret_cast<fftw_complex *>(d_out), FFTW_BACKWARD, FFTW_ESTIMATE);
}

DFT::~DFT()
{
	fftw_free(d_in);
	fftw_free(d_out);
	fftw_destroy_plan(d_plan_fwd);
	fftw_destroy_plan(d_plan_bwd);

	delete[] d_dim;
}

void DFT::forward() { fftw_execute(d_plan_fwd); }
void DFT::backward() { fftw_execute(d_plan_bwd); }
std::complex<double> *DFT::in() { return d_in; }
std::complex<double> *DFT::out() { return d_out; }

