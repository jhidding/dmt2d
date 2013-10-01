#pragma once

#include "complex_util.h"
#include "conan/cube.h"
#include "conan/boxconfig.h"
#include "conan/kspace.h"
#include "dft/dft.h"

using Conan::mVector;
using Conan::BoxConfig;

template <int R>
void calc_potential(DFT &dft, complex64 const *d_hat, Conan::Cube<double> &P)
{
	int const 	N 	= BoxConfig::N();
	double const 	L 	= BoxConfig::L();
	size_t const 	size 	= BoxConfig::size();
	Conan::KSpace<R> K(N, L);		// K-space

	std::transform(d_hat, d_hat + size, K.begin(), dft.in(),
		[] (complex64 const &z, mVector<double, R> const &k)
		{ return z / Conan::sqr(k); });

	dft.in()[0] = 0;
	dft.backward();

	std::transform(dft.out(), dft.out() + size, P.begin(),
		[size] (complex64 const &z)
		{ return z.real() / size; });
}

template <int R>
void calc_hessian_component(DFT &dft, Conan::KSpace<R> const &K, 
		complex64 const *d_hat, double *H, int i, int j)
{
	size_t const 	size 	= BoxConfig::size();

	std::transform(d_hat, d_hat + size, K.begin(), dft.in(),
		[i, j] (complex64 const &z, mVector<double, R> const &k)
		{ return z * (k[i] * k[j]) / Conan::sqr(k); });

	dft.in()[0] = 0;
	dft.backward();

	std::transform(dft.out(), dft.out() + size, H,
		[size] (complex64 const &z)
		{ return z.real() / size; });
}

template <int R>
void calc_hessian(DFT &dft, complex64 const *d_hat, Conan::Cube<double> *H);

template <>
void calc_hessian<2>(DFT &dft, complex64 const *d_hat, Conan::Cube<double> *H)
{
	int const 	N 	= BoxConfig::N();
	double const 	L 	= BoxConfig::L();
	Conan::KSpace<2> K(N, L);		// K-space
	calc_hessian_component(dft, K, d_hat, H[0].data(), 0, 0);
	calc_hessian_component(dft, K, d_hat, H[1].data(), 0, 1);
	calc_hessian_component(dft, K, d_hat, H[2].data(), 1, 1);
}

template <>
void calc_hessian<3>(DFT &dft, complex64 const *d_hat, Conan::Cube<double> *H)
{
	int const 	N 	= BoxConfig::N();
	double const 	L 	= BoxConfig::L();
	Conan::KSpace<3> K(N, L);		// K-space
	calc_hessian_component(dft, K, d_hat, H[0].data(), 0, 0);
	calc_hessian_component(dft, K, d_hat, H[1].data(), 0, 1);
	calc_hessian_component(dft, K, d_hat, H[2].data(), 0, 2);
	calc_hessian_component(dft, K, d_hat, H[3].data(), 1, 1);
	calc_hessian_component(dft, K, d_hat, H[4].data(), 1, 2);
	calc_hessian_component(dft, K, d_hat, H[5].data(), 2, 2);
}

template <int R>
void calc_displacement(DFT &dft, complex64 const *d_hat, Conan::Cube<double> *psi)
{
	int const 	N 	= BoxConfig::N();
	double const 	L 	= BoxConfig::L();
	size_t const 	size 	= BoxConfig::size();
	Conan::KSpace<R> K(N, L);		// K-space

	for (unsigned i = 0; i < R; ++i)
	{
		std::transform(d_hat, d_hat + size, K.begin(), dft.in(),
			[i] (complex64 const &z, mVector<double, R> const &k)
			{ return z * math_i * sin(k[i]); });

		dft.backward();

		std::transform(dft.out(), dft.out() + size, psi[i].begin(),
			[size] (complex64 const &z)
			{ return z.real() / size; });
	}
}

