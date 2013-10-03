/* fourier.h
 *
 * contains methods for fourier space filtering
 */

#pragma once
#include <functional>
#include <cmath>
#include "fft.hh"
#include "mvector.hh"
#include "mdrange.hh"
#include "map.hh"

namespace Fourier
{
	typedef std::complex<double> complex64;
	extern complex64 math_i;

	inline std::function<double (size_t)> wave_number(size_t N, double L)
	{
		return [N, L] (size_t i)
		{
			long a(i);
			return (i > N/2 ? a - long(N) : a) * (2 * M_PI / L);
		};
	}

	inline std::function<double (std::complex<double> const &)> real_part(double size)
	{
		return [size] (std::complex<double> const &z)
		{
			return z.real() / size;
		};
	}

	template <unsigned R>
	using KSpace = System::Map<System::MdRange<R>, System::mVector<double, R>>;

	template <unsigned R>
	KSpace<R> kspace(unsigned N, double L)
	{
		typename KSpace<R>::arg_type shape(N);
		System::MdRange<R> X(shape);

		return KSpace<R>(X, [N, L] (typename KSpace<R>::arg_type const &x)
		{ 
			typename KSpace<R>::value_type k; 
			for (unsigned i = 0; i < R; ++i) 
				k[i] = (unsigned(x[i]) > N/2 ? x[i] - long(N) : x[i]) * (2 * M_PI / L);
			return k; 
		});
	}

	template <unsigned R>
	class Fourier
	{
		public:
			typedef System::mVector<double, R> Vector;
			typedef std::function<double (double)> RealFunction;

			typedef std::function<complex64 (Vector const &)> FilterBase;
			class Filter: public FilterBase
			{
				public:
					Filter(FilterBase const &f):
						FilterBase(f) {}

					Filter operator*(Filter const &o) const
					{
						Filter p(*this);
						return Filter([p, o] (Fourier<R>::Vector const &K)
						{
							return p(K) * o(K);
						});
					}
			};

			static inline Filter scale(double t)
			{
				return Filter([t] (Vector const &K)
				{
					double v = 1.0;
					for (auto k : K)
					{
						v *= exp(t * (cos(k) - 1));
					}
					return complex64(v);
				});
			}

			static inline Filter power_spectrum(RealFunction const &P)
			{
				return Filter([P] (Vector const &K)
				{
					return complex64(sqrt(P(K.norm())));
				});
			}

			static inline Filter potential()
			{
				return Filter([] (Vector const &K)
				{
					return complex64(-1. / K.sqr());
				});
			}

			static inline Filter derivative(unsigned i)
			{
				return Filter([i] (Vector const &K)
				{
					return math_i * sin(K[i]);
				});
			}

			static inline std::function<complex64 (complex64, Vector const &)>
			filter(Filter const &f)
			{
				return [f] (complex64 z, Vector const &K)
				{
					return z * f(K);
				};
			}
	};
}

// vim:ts=4:sw=4:tw=80
