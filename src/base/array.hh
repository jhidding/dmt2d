/* array.h
 *
 * Arrays are vectors with a shape tuple
 */

#pragma once
#include "misc.hh"
#include "mtypeid.hh"
#include "header.hh"
#include <iostream>
#include <vector>
#include <iterator>

namespace System
{
	// Ref {{{1
	/*!
	 * reference counted pointer to a range-type. maybe up for a
	 * name change, but should remain short.
	 */
	template <typename T>
	class Ref
	{
		std::shared_ptr<T> _data;

		public:
			typedef typename T::value_type value_type;
			typedef typename T::iterator iterator;
			typedef typename T::const_iterator const_iterator;

			Ref() {}

			Ref(T *data_): 
				_data(data_) 
			{}

			T &operator*() { return *_data; }
			T *operator->() { return _data.get(); }
			T *get() { return _data.get(); }

			size_t size() const { return _data->size(); }
			value_type &operator[](size_t idx) { return (*_data)[idx]; }
			value_type const &operator[](size_t idx) const { return (*_data)[idx]; }
			iterator begin() const { return _data->begin(); }
			iterator end() const { return _data->end(); }	
			const_iterator cbegin() const { return _data->cbegin(); }
			const_iterator cend() const { return _data->cend(); }	
	};
	// }}}1

	template <typename T>
	class Array: public Ref<std::vector<T>>
	{
		public:
			using Ref<std::vector<T>>::get;

			Array(size_t n, T value = T()):
				Ref<std::vector<T>>(new std::vector<T>(n, value))
			{}

			Array(std::istream &in):
				Ref<std::vector<T>>(new std::vector<T>)
			{
				read_block(in, *get());
			}

			void to_file(std::ostream &out)
			{
				write_block(out, *get());
			}
	};

	template <typename T>
	void save_to_file(std::ostream &fo, Array<T> data, std::string const &name)
	{
		Header S;
		S["name"] = name;
		S["dtype"] = TypeRegister::name<T>();
		S.to_file(fo); data.to_file(fo);
	}

	template <typename T>
	Array<T> load_from_file(std::istream &fi, std::string const &name)
	{
		auto pos = fi.tellg();

		while (fi.good())
		{
			Header S(fi);
			if (S["name"] != name)
			{
				skip_block(fi);
			}
			else
			{
				Array<T> data(fi);
				fi.seekg(pos, std::ios::beg);
				return data;
			}
		}

		throw "Couldn't find record " + name + " in file.";
	}
}

// vim:ts=4:sw=4:tw=80:fdm=marker
