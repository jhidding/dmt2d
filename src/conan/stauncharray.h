/*
 * StaunchArray is called that way, because it is an array that should not be
 * moved and cannot be copied.  It does not use refference counting or other
 * 'bells and whistles'. Just don't let them get out of scope.  The data can be
 * saved or loaded from fortran-like data blocks, used in many scientific
 * applications.
 */

#pragma once

#include <iostream>
#include <fstream>

namespace Conan
{

template <typename T>
class StaunchArray
{
	T		*d;
	unsigned	M;
	unsigned	byte_size;

	public:
		typedef T 	*iterator;
		typedef T const *const_iterator;
		typedef T 	Type;

		StaunchArray()
			: d(NULL)
		{}

		StaunchArray(unsigned M_)
			: M(M_), byte_size(M_ * sizeof(T))
		{
			d = new T[M];
		}

		void destroy()
		{
			if (d != NULL)
				delete[] d;
		}

		~StaunchArray() { destroy(); } 

		StaunchArray(std::istream &f)
		{
			unsigned bs1, bs2;	
			f.read(reinterpret_cast<char *>(&bs1), sizeof(unsigned));
			if (f.fail())
			{
				throw "Could not read from file";
			}
			byte_size = bs1;
			M = byte_size / sizeof(T);
			d = new T[M];

			f.read(reinterpret_cast<char *>(d), byte_size);
			f.read(reinterpret_cast<char *>(&bs2), sizeof(unsigned));

			if (bs1 != bs2)
				throw "The block seems to be corrupted";
		}

		void resize(size_t M_)
		{
			destroy();
			M = M_;
			byte_size = M_ * sizeof(T);
			d = new T[M];
		}

		size_t size() const { return M; } 


		void to_file(std::ostream &f)
		{
			f.write(reinterpret_cast<char *>(&byte_size), sizeof(unsigned));
			f.write(reinterpret_cast<char *>(d), byte_size);
			f.write(reinterpret_cast<char *>(&byte_size), sizeof(unsigned));
		}

		void to_raw_file(std::ostream &f)
		{
			f.write(reinterpret_cast<char *>(d), byte_size);
		}

		void from_file(std::istream &f)
		{
			if (d == NULL)
				throw "You need to allocate memory first.";

			unsigned bs1, bs2;	
			f.read(reinterpret_cast<char *>(&bs1), sizeof(unsigned));

			if (bs1 != byte_size)
				throw "Block size does not agree with my expectations";

			f.read(reinterpret_cast<char *>(d), byte_size);
			f.read(reinterpret_cast<char *>(&bs2), sizeof(unsigned));

			if (bs1 != bs2)
				throw "The block seems to be corrupted";
		}

		Type &operator[](size_t idx) { return d[idx]; } 
		Type const &operator[](size_t idx) const { return d[idx]; } 

		iterator data() { return d; } 
		iterator begin() { return d; } 
		iterator end() { return d + M; } 

		const_iterator data() const { return d; } 
		const_iterator begin() const { return d; } 
		const_iterator end() const { return d + M; } 

	private:
		StaunchArray(StaunchArray const &); // NI
};

} // namespace Conan

