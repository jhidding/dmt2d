/* mvector.h
 *
 * mini vector
 */

#pragma once

#include <algorithm>
#include <functional>
#include <numeric>
#include <initializer_list>
#include "misc.hh"

namespace System
{
	template <typename T, unsigned R>
	class mVector
	{
		T	p[R];

		public:
			typedef T value_type;

			mVector() {}
			mVector(T v)
			{
				std::fill_n(p, R, v);
			}

			template <typename U>
			mVector(mVector<U, R> const &o)
			{
				std::copy(o.begin(), o.end(), begin());
			}

			mVector(std::initializer_list<T> lst)
			{
				std::copy(lst.begin(), lst.end(), p);
			}

			template <typename U>
			mVector<U, R> as() const
			{
				mVector<U, R> v;
				std::copy(begin(), end(), v.begin());
				return v;
			}

			T *begin() { return p; }
			T *end() { return p + R; }
			T const *begin() const { return p; }
			T const *end() const { return p + R; }
			T &operator[](unsigned i) { return p[i]; }
			T const &operator[](unsigned i) const { return p[i]; }

			template <typename Fun>
			mVector operate(Fun f) const
			{ 
				mVector a;
				std::transform(begin(), end(), a.begin(), f);
				return a;
			}

			template <typename Fun>
			mVector operate(Fun f, mVector const &o) const
			{ 
				mVector a;
				std::transform(begin(), end(), o.begin(), a.begin(), f);
				return a;
			}

			mVector operator+=(mVector const &o)
			{ std::transform(begin(), end(), o.begin(), begin(), std::plus<T>()); return *this; }
			mVector operator-=(mVector const &o)
			{ std::transform(begin(), end(), o.begin(), begin(), std::minus<T>()); return *this; }
			mVector operator*=(mVector const &o)
			{ std::transform(begin(), end(), o.begin(), begin(), std::multiplies<T>()); return *this; }

			template <typename U>
			mVector operator/=(U o)
			{ 
				std::transform(begin(), end(), begin(),
						std::bind2nd(std::divides<T>(), o)); 
				return *this; 
			}
			template <typename U>
			mVector operator*=(U o)
			{ 
				std::transform(begin(), end(), begin(),
						std::bind2nd(std::multiplies<T>(), o)); 
				return *this; 
			}
			template <typename U>
			mVector operator%=(U o)
			{ 
				std::transform(begin(), end(), begin(), 
						[o] (T a) { return System::modulus(a, o); }); 
				return *this; 
			}

			mVector operator+(mVector const &o) const
			{ return operate(std::plus<T>(), o); }
			mVector operator-(mVector const &o) const
			{ return operate(std::minus<T>(), o); }
			mVector operator-() const
			{ return operate(std::negate<T>()); }
			mVector operator*(mVector const &o) const
			{ return operate(std::multiplies<T>(), o); }

			template <typename U>
			mVector operator/(U o) const
			{ return operate(std::bind2nd(std::divides<T>(), o)); }
			template <typename U>
			mVector operator*(U o) const
			{ return operate(std::bind2nd(std::multiplies<T>(), o)); }
			template <typename U>
			mVector operator%(U o) const
			{ return operate([o] (T a) { return System::modulus(a, o); }); }

			// maths
			mVector abs() const
			{
				mVector a; 
				std::transform(begin(), end(), a.begin(), ::fabs);
				return a;
			}

			T sqr() const
			{ 
				T v = 0; 
				for (auto x : *this)
					v += x*x;
				return v;
			}

			T norm() const
			{ 
				return sqrt(sqr());
			}

			T max() const
			{
				return *std::max_element(begin(), end());
			}

			T sum() const
			{
				return std::accumulate(begin(), end(), T(0));
			}

			T inf_norm() const
			{
				return abs().max();
			}

			T dot(mVector const &o) const
			{
				return std::inner_product(begin(), end(), o.begin(), T(0));
			}
	};

	template <typename T, unsigned R>
	std::ostream &operator<<(std::ostream &out, mVector<T, R> const &v)
	{
		std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
		return out;
	}

}

// vim:ts=4:sw=4:tw=80
