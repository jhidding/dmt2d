#pragma once

#include <iterator>
#include "mvector.h"

/*!
 * \ingroup spaces
 * A Space iterating over an array of integer vectors in the order
 * of array indexing.
 */

namespace Conan {

template <int rank>
class NdSpace
{
	mVector<int, rank>	shape;
	size_t				n;

	public:
		typedef mVector<int, rank> value_type;

		class const_iterator;

		NdSpace(value_type const &shape_)
			: shape(shape_)
		{
			n = std::accumulate(
				shape.begin(), shape.end(), size_t(1),
				std::multiplies<size_t>());
		}

		value_type operator[](size_t i) const
		{
			value_type q;
			//for (int k = rank-1; k >= 0; --k)
			for (int k = 0; k < rank; ++k)
			{
				q[k] = i % shape[k];
				i /= shape[k];
			}
			return q;
		}

		const_iterator begin() const
		{
			return const_iterator(shape, value_type(int(0)), 0);
		}

		const_iterator end() const
		{
			return const_iterator(shape, shape, n);
		}
};

template <int rank>
class NdSpace<rank>::const_iterator: public std::iterator<std::forward_iterator_tag, const value_type>
{
	mVector<int, rank> shape, value;
	size_t			n;

	public:
		const_iterator(value_type const &shape_, value_type const &value_, size_t n_)
			: shape(shape_), value(value_), n(n_)
		{
		}

		const_iterator &operator++()
		{
			++value[0];
			++n;
			cascade();

			//if (value[rank-1] >= shape[rank-1])
			//{
			//	value[rank-1] = 0;
			//}

			return *this;
		}

		const_iterator &operator+=(size_t a)
		{
			value[0] += a;
			n += a;
			cascade();
			return *this;
		}

		const_iterator operator+(size_t a)
		{
			return const_iterator(*this) += a;
		}

		void cascade()
		{
			for (int i = 0; i < rank - 1; ++i)
			{
				if (value[i] >= shape[i])
				{
					value[i+1] += value[i] / shape[i];
					value[i] %= shape[i];
				}
				else
				{
					break;
				}
			}
		}

		value_type const &operator*() const
			{ return value; }

		bool operator==(const_iterator const &other) const
			{ return n == other.n; }

		bool operator!=(const_iterator const &other) const
			{ return n != other.n; }
};

}

