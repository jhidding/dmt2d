#pragma once

#include <iostream>
#include <iterator>
#include <algorithm>

#include "stauncharray.h"
#include "boxconfig.h"
#include "shape_str.h"
#include "header.h"
#include "mtypeid.h"
#include "discard_block.h"
#include "convert.h"

namespace Conan
{
	template <typename T>
	class Cubelet: public StaunchArray<T>
	{
		public:
			unsigned	rank, bits;
			unsigned	N;

			Cubelet(unsigned N_):
				StaunchArray<T>(1L << (BoxConfig::rank() * log_binary(N_))),
				rank(BoxConfig::rank()), bits(log_binary(N_)),
				N(1 << bits)
			{}
	};

	template <typename T>
	class Cube: public StaunchArray<T>
	{
		std::string d_name;

		public:
			typedef T type;

			Cube():
				StaunchArray<T>(BoxConfig::size())
			{}

			Cube(std::string const &name_):
				StaunchArray<T>(BoxConfig::size()),
				d_name(name_)
			{}

			Cube(std::istream &in);
			void rename(std::string const &name_);
			std::string name() const;
			void to_file(std::ostream &out);
			void from_file(std::istream &in);
			void convert_from_file(std::istream &in);
			bool search_in_file(std::istream &fi, std::string const &name);
	};

	template <typename T>
	class CubeSlice: public std::iterator<typename T::type, 
					std::forward_iterator_tag>
	{
		typename Cube<T>::iterator p;
		unsigned i;

		public:
			typedef typename T::type type;

			CubeSlice(typename Cube<T>::iterator p_, unsigned i_):
				p(p_), i(i_) {}

			CubeSlice &operator++()
			{ ++p; return *this; }

			type &operator*()
			{ return (*p)[i]; }

			bool operator!=(CubeSlice const &o) const
			{ return p != p.o; }
	};

	template <typename T>
	bool Cube<T>::search_in_file(std::istream &fi, std::string const &name)
	{
		while (true)
		{
			try 
			{
				std::streampos pos = fi.tellg();
				Header H(fi);

				if (H["name"] == name)
				{
					fi.seekg(pos);
					convert_from_file(fi);
					return true;
				}
				else
				{
					discard_block(fi);
				}
			}
			catch (char const *c)
			{
				break;
			}
		}

		return false;
	}

	template <typename T>
	std::string Cube<T>::name() const
	{
		return d_name;
	}

	template <typename T>
	Cube<T>::Cube(std::istream &in):
		StaunchArray<T>(BoxConfig::size())
	{
		from_file(in);
	}

	template <typename T>
	void Cube<T>::rename(std::string const &name_)
	{
		d_name = name_;
	}

	template <typename T>
	void Cube<T>::to_file(std::ostream &out)
	{
		Header H;
		H["type"]  = mTypeId<T>();
		H["shape"] = shape_str(BoxConfig::N(), BoxConfig::rank());
		H["N"] = to_string(BoxConfig::N());
		H["L"] = to_string(BoxConfig::L());
		H["rank"] = to_string(BoxConfig::rank());
		H["name"] = d_name;

		H.to_file(out);
		StaunchArray<T>::to_file(out);
	}

	template <typename T>
	void Cube<T>::from_file(std::istream &in)
	{
		Header H(in);
		if (H["type"] != mTypeId<T>())
			throw "Type mismatch in reading from file";

		if (from_string<unsigned>(H["rank"]) != BoxConfig::rank())
			throw "Wrong dimension, please try the eleventh. (the errors you get these days...)";

		if (from_string<unsigned>(H["N"]) != BoxConfig::N())
			throw "Wrong box size";

		d_name = H["name"];

		StaunchArray<T>::from_file(in);
	}

	template <typename T>
	void Cube<T>::convert_from_file(std::istream &in)
	{
		Header H(in);
	
		if (from_string<unsigned>(H["rank"]) != BoxConfig::rank())
			throw "Wrong dimension, please try the eleventh. (the errors you get these days...)";

		if (from_string<unsigned>(H["N"]) != BoxConfig::N())
			throw "Wrong box size";

		d_name = H["name"];
		if (H["type"] == mTypeId<T>())
		{
			std::cerr << "Reading block in native format: " << mTypeId<T>() << std::endl;
			StaunchArray<T>::from_file(in);
			return;
		}

		std::cerr << "Warning: attempting to convert data from " << H["type"] << " to " 
			  << mTypeId<T>() << "." << std::endl;

		if (H["type"] == mTypeId<int>())
		{
			StaunchArray<int> temp(in);
			std::copy(temp.begin(), temp.end(), StaunchArray<T>::begin());
		}

		if (H["type"] == mTypeId<float>())
		{
			StaunchArray<float> temp(in);
			std::copy(temp.begin(), temp.end(), StaunchArray<T>::begin());
		}

		if (H["type"] == mTypeId<double>())
		{
			StaunchArray<double> temp(in);
			std::copy(temp.begin(), temp.end(), StaunchArray<T>::begin());
		}
	}

}

