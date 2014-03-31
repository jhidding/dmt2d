#pragma once

#include <vector>
#include "../base/misc.hh"
#include "../base/array.hh"
#include "../base/fourier.hh"
#include "../base/fft.hh"

namespace DMT {
	using System::Array;
	using System::copy;
	using System::BoxConfig;
	using System::ptr;
	using System::complex64;

	template <unsigned R>
	class HessianBase: public Array<Array<double>>
	{
		// std::vector<Array<double>> data;

		public:
			HessianBase(ptr<BoxConfig<R>> box, Array<double> A);

			// Array<double> operator[](unsigned i) const { return data[i]; }
	};

	template <unsigned R>
	HessianBase<R>::HessianBase(ptr<BoxConfig<R>> box, Array<double> A):
		Array<Array<double>>(0,0)
	{
		for (unsigned i = 0; i < (R * (R + 1))/2; ++i)
			get()->push_back(Array<double>(box->size()));

		Fourier::Transform dft(box->shape());
		Array<complex64> A_f(box->size());
		copy(A, dft.in);
		dft.forward();
		copy(dft.out, A_f);

		auto K = Fourier::kspace<R>(box->N(), box->L());
		unsigned o = 0;
		double s = box->size() / box->scale2();	
		for (unsigned i = 0; i < R; ++i) {
			for (unsigned j = 0; j <= i; ++j)
		{
			//push_back(Array<double>(box->size()));
			auto F = Fourier::Fourier<R>::derivative(i) 
			       * Fourier::Fourier<R>::derivative(j);
			transform(A_f, K, dft.in, Fourier::Fourier<R>::filter(F));
			dft.backward();
			transform(dft.out, (*this)[o], Fourier::real_part(box->size()));
			++o;
		} }
	}

	template <unsigned R>
	class Hessian;

	template <>
	class Hessian<2>: public HessianBase<2>
	{
		public:
			using HessianBase::HessianBase;
	};
}

