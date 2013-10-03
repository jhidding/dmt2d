/* fft.h
 *
 * compute the Fast Fourier transform using FFTW
 */

#pragma once

#include <complex>
#include <vector>
#include <algorithm>
#include <fftw3.h>
#include <functional>
#include <iterator>
#include <memory>

namespace Fourier
{
	template <typename T>
	class FFT_allocator: public std::allocator<T>
	{
		public:
			typedef T		value_type;
			typedef T *	 	pointer;
			typedef T &	 	reference;
			typedef T const * 	const_pointer;
			typedef T const &	const_reference;
			typedef size_t 		size_type;
			typedef ptrdiff_t 	difference_type;

			pointer allocate(size_t n, std::allocator<void>::const_pointer hint = 0)
			{
				if (hint != 0)
					fftw_free(hint);

				return reinterpret_cast<T *>(fftw_malloc(n * sizeof(T)));
			}

			void deallocate(pointer p, size_t n)
			{
				fftw_free(p);
			}
	};

	class Transform
	{	
		size_t			size;
		fftw_plan		d_plan_fwd, d_plan_bwd;

		public:
			std::vector<std::complex<double>, FFT_allocator<std::complex<double>>> in, out;

			Transform(std::vector<int> const &);
			~Transform();
			void forward();
			void backward();
	};
}

