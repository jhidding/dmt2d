#include "ic.hh"
#include "../base/fourier.hh"
#include "../base/access.hh"

using namespace System;

template <typename A>
typename A::value_type mean(A const &a)
{
	return std::accumulate(a.begin(), a.end(), typename A::value_type(0)) / a.size();
}

template <unsigned R>
Array<double> _generate_random_field(Header const &C)
{
	unsigned    mbits = C.get<unsigned>("mbits");
	double	        L = C.get<double>("size");
	double 	    slope = C.get<double>("power-slope");
	bool       smooth = C.get<bool>("smooth");
	double      sigma = C.get<double>("scale");
	unsigned     seed = C.get<unsigned>("seed");

	size_t N = size_t(1) << mbits;

	mVector<int, R> shape(N);
	size_t size = product(shape);
	Fourier::Transform fft(std::vector<int>(R, N));

	Array<double> dens(size);
	generate(dens, Gaussian_white_noise(seed));
	copy(dens, fft.in);

	auto P = Fourier::Fourier<R>::power_spectrum(
		[slope] (double k) { return pow(k, slope); });
	auto S = Fourier::Fourier<R>::scale(sigma * N/L);
	auto K = Fourier::kspace<R>(N, N);

	fft.forward();
	transform(fft.out, K, fft.in, Fourier::Fourier<R>::filter(P * S));
	fft.in[0] = 0;
	fft.backward();
	copy(dens, fft.in);
	transform(fft.out, dens, Fourier::real_part(size));

	double var = mean(map(dens, [] (double a) { return a*a; }));

	fft.forward();
	transform(fft.out, K, fft.in, 
			Fourier::Fourier<R>::filter((smooth ? P * S : P)));

	fft.in[0] = 0;
	fft.backward();
	transform(fft.out, dens, Fourier::real_part(size * sqrt(var)));
	
	return dens;
}

template <unsigned R>
void _compute_potential(Header const &C, Array<double> density)
{
	double L = C.get<double>("size");
	unsigned mbits = C.get<unsigned>("mbits");
	unsigned N = 1 << mbits;
	size_t size = 1U << (mbits * R);

	Fourier::Transform fft(std::vector<int>(R, N));
	auto K = Fourier::kspace<R>(N, L);
	auto F = Fourier::Fourier<R>::potential();
	copy(density, fft.in);
	fft.forward();
	transform(fft.out, K, fft.in, Fourier::Fourier<R>::filter(F));
	fft.in[0] = 0;
	fft.backward();
	transform(fft.out, density, Fourier::real_part(size));
}

template <unsigned R>
Array<mVector<double, R>> _compute_displacement(Header const &C, Array<double> potential)
{
	double L = C.get<double>("size");
	unsigned mbits = C.get<unsigned>("mbits");
	unsigned N = 1 << mbits;
	size_t size = 1U << (mbits * R);

	Array<Fourier::complex64> phi_f(size);
	Fourier::Transform fft(std::vector<int>(R, N));
	auto K = Fourier::kspace<R>(N, N);
	copy(potential, fft.in);
	fft.forward();
	copy(fft.out, phi_f);

	Array<mVector<double, R>> psi(size);
	for (unsigned k = 0; k < R; ++k)
	{
		auto F = Fourier::Fourier<R>::derivative(k);
		auto psi_k = access(psi, [k] (mVector<double, R> &x) -> double&
			{ return x[k]; });

		transform(phi_f, K, fft.in, Fourier::Fourier<R>::filter(F));
		fft.backward();
		transform(fft.out, psi_k, Fourier::real_part(size / L * N));
	}
	return psi;
}

Array<double> Conan::generate_random_field(Header const &C)
{
	unsigned      dim = C.get<unsigned>("dim");

	switch (dim)
	{
		case 2: return _generate_random_field<2>(C);
		case 3: return _generate_random_field<3>(C);
	}

	throw "only 2 and 3 dimensions supported.";
}

void Conan::compute_potential(Header const &C, Array<double> data)
{
	unsigned      dim = C.get<unsigned>("dim");

	switch (dim)
	{
		case 2: _compute_potential<2>(C, data); return;
		case 3: _compute_potential<3>(C, data); return;
	}

	throw "only 2 and 3 dimensions supported.";
}

void Conan::compute_displacement(Header const &C, Array<double> data, std::ostream &fo)
{
	unsigned dim = C.get<unsigned>("dim");

	switch (dim)
	{
		case 2: save_to_file(fo, _compute_displacement<2>(C, data), "displacement");
			return;
		case 3: save_to_file(fo, _compute_displacement<3>(C, data), "displacement");
			return;
	}
	throw "only 2 and 3 dimensions supported.";
}

