#pragma once
#include "mvector.hh"
#include "cvector.hh"

namespace System
{
	/*!
	 * struct containing data members giving the box paramaters of
	 * the simulation. These describe the dimension, simulation size and physical size
	 * of the simulation
	 */
	template <unsigned R>
	class BoxConfig
	{
		cVector<R>	d_box;
		std::vector<int> d_shape;
		unsigned 	d_N;
		unsigned	d_bits;
		size_t		d_size;
		double		d_L;
		double		d_scale, d_scale2;

		public:
			BoxConfig(unsigned bits_, float L_);	///< set the values, calculates bits and size
			unsigned N() const { return d_N; } 		///< The number of pixels on each side
			unsigned bits() const { return d_bits; } 	///< The number of bits = 2log(N)
			size_t size() const { return d_size; }	///< The size of the box = N^rank
			double L() const { return d_L; }		///< Length of the box in units of Mpc h^-1
			double scale() const { return d_scale; }
			double scale2() const { return d_scale2; }
			std::vector<int> const &shape() const { return d_shape; }

			cVector<R> const &box() const { return d_box; }
			/*typename cVector<R>::C m2c(mVector<int, R> const &v) const
			{ return d_box.loc(v); }
			std::vector<typename cVector<R>::C> const &sq() const
			{ return d_box.sq; }*/
	};

	/// calculate the binary exponent i * 2^shift
	inline size_t pow_binary(size_t i, unsigned shift)
	{
		return i << shift;
	}

	template <unsigned R>
	BoxConfig<R>::BoxConfig(unsigned bits_, float L_):
		d_box(bits_), d_N(pow_binary(1, bits_)), d_bits(bits_),
		d_size(pow_binary(1, R * d_bits)),
		d_L(L_), d_scale(d_L/d_N), 
		d_scale2(d_scale * d_scale)
	{
		for (unsigned i = 0; i < R; ++i) d_shape.push_back(d_N);
	}
}

