#pragma once

#include <functional>
#include <iterator>

namespace Conan {

template <typename In, typename Out>
class Map
{
	typedef std::function<Out (typename In::value_type)> F;
	F 	 fn;
	In const &in;

	public:
		typedef Out value_type;

		class const_iterator: public std::iterator<std::forward_iterator_tag, const value_type>
		{
			F 				fn;
			typename In::const_iterator	i;

			public:
				const_iterator(F const &fn_, typename In::const_iterator const &i_): 
					fn(fn_), i(i_) {}

				const_iterator &operator++()
					{ ++i; return *this; }

				const_iterator &operator--()
					{ --i; return *this; }

				const_iterator &operator+=(size_t a)
					{ i += a; return *this; }

				bool operator!=(const_iterator const &other) const
					{ return i != other.i; }

				bool operator==(const_iterator const &other) const
					{ return i == other.i; }

				value_type operator*() const { return fn(*i); }
		};

		Map(F const &fn_, In const &in_):
			fn(fn_), in(in_) {}
		
		value_type operator[](size_t i)
			{ return fn(in[i]); }

		const_iterator begin() const 
			{ return const_iterator(fn, in.begin()); }
		const_iterator end() const 
			{ return const_iterator(fn, in.end()); }
};

}
