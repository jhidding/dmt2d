#pragma once

/*!
 * \class mVector
 * Simplifies a lot of vector calculations. the concept is
 * that the vector contains no other data then its values
 * the meta data is all handled by the c++ template system
 * this makes the mVector very easy to use.
 *
 * A few auxiliary functions are defined (and behave as 
 * expected from their names):
 * 	map, abs, dot
 *
 * also the ostream insertion operator is overloaded, so that
 * the vector can be printed on screen using tab as delimeter.
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <functional>
#include <iterator>
#include <cmath>

#include "modulus.h"

namespace Conan
{

template <typename T, typename U>
T	cast_function(U const &value)
{
	return static_cast<T>(value);
}

inline double _round_up(double a)
{ return (a < 0.0 ? int(a) : int(a) + 1); }

inline double _round_down(double a)
{ return (a < 0.0 ? int(a) - 1: int(a)); }

/*!
 * obsolete: should use bind2nd, but C++ can't resolve which modulus
 * to use as the function is overloaded, should solve this.
 */
template <typename T>
class modder
{
	T a;
	public:
		modder(T a_): a(a_) {}
		T operator()(T b) { return modulus(b, a); }
};

template <typename T, int rank>
class mVector
{
	protected:
		T	d[rank];

	public:
		enum { dimension = rank };
		typedef T type;
		typedef T* iterator;

		mVector(T def = T(0))
			{ std::fill_n(d, rank, def); }

//		mVector(T *v)
//			{ std::copy(v, v + rank, d); }

		template <typename U>
		mVector(mVector<U, rank> const &other)
		{
			std::transform(other.data(), other.data() + rank, d, cast_function<T,U>);
		}
		// casts -----------------------------------------
		template <typename U>
		operator mVector<U, rank>()
		{
			return mVector<U, rank>(*this);
		}

		// vector addition -------------------------------
		template <typename U>
		mVector &operator+=(mVector<U, rank> const &other)
		{
			std::transform(d, d + rank, other.data(), d, std::plus<T>());
			return *this;
		}

		template <typename U>
		mVector operator+(mVector<U, rank> const &other) const
		{
			mVector<T, rank> x(*this);
			x += other;
			//for (int i = 0; i < rank; ++i) x[i] = d[i] + other[i];
			return x;
			//return mVector(*this) += other;
		}

		mVector operator+(T a) const
		{
			return *this + mVector(a);
		}

		// vector substraction -------------------------
		template <typename U>
		mVector &operator-=(mVector<U, rank> const &other)
		{	
			std::transform(d, d + rank, other.data(), d, std::minus<T>());
			return *this;
		}

		template <typename U>
		mVector operator-(mVector<U, rank> const &other) const
			{ return mVector(*this) -= other; }

		mVector operator-(T a) const
		{
			return *this - mVector(a);
		}

		mVector operator-() const
		{
			return (*this) * (-1);
		}

		// scalar multiplication -------------------------------
		mVector &operator*=(T scalar)
		{
			std::transform(d, d + rank, d, 
				std::bind1st(std::multiplies<T>(), scalar));
			return *this;
		}

		template <typename U>
		mVector operator*(U scalar) const
			{ return mVector(*this) *= scalar; }

		// scalar division -------------------------------
		mVector &operator/=(T scalar)
		{
			std::transform(d, d + rank, d, 
				std::bind2nd(std::divides<T>(), scalar));
			return *this;
		}

		template <typename U>
		mVector operator/(U scalar) const
			{ return mVector(*this) /= scalar; }

		mVector inv() const
		{
			mVector u;
			std::transform(d, d + rank, u,
				std::bind1st(std::divides<T>(), 1.0));
			return u;
		}

		// modulus of single scalar, or piecewise vector --------
		mVector &operator%=(mVector const &other)
		{
			std::transform(d, d+rank, other.d, d, Conan::modulus);
			return *this;
		}

		mVector &operator%=(T scalar)
		{
			std::transform(d, d + rank, d, modder<T>(scalar));
			return *this;
		}

		mVector operator%(mVector<T, rank> const &other) const
			{ return mVector(*this) %= other; }

		mVector operator%(T scalar) const
			{ return mVector(*this) %= scalar; }

		/*
		 * rather explicit version of the modulus function, probably made
		 * when the above members didn't function as expected for some unclarified
		 * reason, this version is depricated
		 */
		//mVector &modulus(T scalar)
		//{
		//	d[0] -= scalar * ::floor(d[0] / scalar);
		//	d[1] -= scalar * ::floor(d[1] / scalar);
		//	d[2] -= scalar * ::floor(d[2] / scalar);
		//	return *this;
		//}

		mVector floor() const
		{
			mVector v;
			std::transform(begin(), end(), v.begin(), _round_down);
			return v;
		}

		mVector ceil() const
		{
			mVector v;
			std::transform(begin(), end(), v.begin(), _round_up);
			return v;
		}

		/*!
		 * reduce vector to one of a rank one less than current, taking
		 * out one value at position k. not very efficient, also don't know
		 * where I'm using this.
		 */
		mVector<T, rank-1> reduce(size_t k) const
		{
			mVector<T, rank - 1> r;
			size_t j = 0;
			for (size_t i = 0; i < rank; i++)
			{
				r[j] = d[i];
				if (i != k) ++j;
			}

			return r;
		}

		// accessor functions ------------------------------------
		T &operator[](int i) { return *(d + i); }
		T const &operator[](int i) const { return *(d + i); }
		T *data() { return d; }
		T const *data() const { return d; }
		T *begin() { return d; }
		T *end() { return d + rank; }
		T const *begin() const { return d; }
		T const *end() const { return d + rank; }
};

template <typename T, typename U, int rank>
T dot(mVector<T, rank> const &a, mVector<U, rank> const &b)
{
	return std::inner_product(a.begin(), a.end(), b.begin(), T(0));
}

template <typename T, typename U, int rank>
T trace(mVector<T, rank> const &a, mVector<U, rank> const &b)
{
	return std::inner_product(a.begin(), a.end(), b.begin(), T(1.0), 
			std::multiplies<T>(), std::multiplies<T>());
}

template <typename T, int rank>
inline double norm(mVector<T, rank> const &v)
{
	return sqrt(std::inner_product(v.begin(), v.end(), v.begin(), T(0)));
}

template <typename T, int rank>
inline double sqr(mVector<T, rank> const &v)
{
	return std::inner_product(v.begin(), v.end(), v.begin(), T(0));
}

template <typename T>
inline mVector<T, 3> cross(mVector<T, 3> const &a, mVector<T, 3> const &b)
{
	mVector<T, 3> c;

	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];

	return c;
}

template <typename T, int rank, typename U>
mVector<T, rank> mmap(U fn, mVector<T, rank> const &v)
{
	mVector<T, rank> a;
	std::transform(v.begin(), v.end(), a.begin(), fn);
	return a;
}

template <typename T, typename V, int rank, typename U>
mVector<T, rank> mmap(U fn, mVector<T, rank> const &v1, mVector<V, rank> const &v2)
{
	mVector<T, rank> a;
	std::transform(v1.begin(), v1.end(), v2.begin(), a.begin(), fn);
	return a;
}

template <typename T, int rank>
std::ostream &operator<<(std::ostream &out, mVector<T, rank> const &v)
{
	std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, "\t"));
	return out;
}

} // namespace Conan

