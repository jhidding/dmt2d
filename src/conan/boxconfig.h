#pragma once

#include "metric.h"

/*!
 * \namespace Conan
 * \brief Conan contains classes that provide a generic interface for the different
 * modules.
 */

namespace Conan
{

/*! 
 * calculate the binary logarithm of N, actually calculates the number of bits
 * one needs to represent N - 1
 */
inline unsigned log_binary(unsigned N)
{
	if (N <= 1) return 0;
	else return log_binary(N >> 1) + 1;
}

/// calculate the binary exponent i * 2^shift
inline size_t pow_binary(size_t i, unsigned shift)
{
	return i << shift;
}

/*!
 * struct containing data members giving the box paramaters of
 * the simulation. These describe the dimension, simulation size and physical size
 * of the simulation
 */
class BoxConfig
{
	static BoxConfig *singleton;

	unsigned 	d_N;
	unsigned	d_bits;
	unsigned	d_rank;
	size_t		d_size;
	float		d_L;
	Metric		d_M;
	float		d_scale, d_scale2;

	private:
		BoxConfig(unsigned rank, unsigned N_, float L_);	///< set the values, calculates bits and size

	public:

		static BoxConfig const &initialize(unsigned rank, unsigned N_, float L_);
		static BoxConfig const &instance();

		static unsigned N() { return singleton->d_N; } 		///< The number of pixels on each side
		static unsigned bits() { return singleton->d_bits; } 	///< The number of bits = 2log(N)
		static unsigned rank() { return singleton->d_rank; } 	///< The number of bits = 2log(N)
		static size_t size() { return singleton->d_size; }	///< The size of the box = N^rank
		static float L() { return singleton->d_L; }		///< Length of the box in units of Mpc h^-1
		static float scale() { return singleton->d_scale; }
		static float scale2() { return singleton->d_scale2; }

		template <typename T>
		static size_t M(T const &v) { return singleton->d_M(v); }
};


} // namespace Conan

