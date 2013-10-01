#pragma once

#include <functional>
#include <iterator>

template <typename T>
class Repeat
{
	typedef std::function<T ()> F;
	F	fn;
	size_t 	N;

	public:
		typedef T value_type;

		class const_iterator: public std::iterator<std::forward_iterator_tag, value_type>
		{
			F 	fn;
			size_t 	i;
			T	value;

			public:
				const_iterator(F const &fn_, size_t i_): 
					fn(fn_), i(i_), value(fn()) {}

				const_iterator &operator++()
					{ ++i; value = fn(); return *this; }

				const_iterator &operator--()
					{ --i; return *this; }

				const_iterator &operator+=(size_t a)
					{ i += a; return *this; }

				bool operator!=(const_iterator const &other) const
					{ return i != other.i; }

				bool operator==(const_iterator const &other) const
					{ return i == other.i; }

				value_type &operator*() const { return value; }
				value_type const *operator->() const { return &value; } 
		};

		Repeat(F const &fn_, size_t N_):
			fn(fn_), N(N_) {}

		const_iterator begin() const 
			{ return const_iterator(fn, 0); }
		const_iterator end() const 
			{ return const_iterator(fn, N); }
};
