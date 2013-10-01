#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

namespace Conan
{

class ByteVector
{
	std::vector<char>		d_data;
	std::vector<char>::iterator	d_ptr;

	template <typename T>
	class convert
	{
		union
		{
			char 	d_raw[sizeof(T)];
			T	d_datum;
		};

		public:
			convert(char *c)
			{
				std::copy(c, c + sizeof(T), d_raw);
			}

			convert(T const &d): d_datum(d) {}

			T const		&datum() const { return d_datum; }
			char const	*raw() const { return d_raw; }
	};

	public:
		ByteVector() {}

		ByteVector(std::istream &fi)
		{
			from_file(fi);
		}

		template <typename T>
		ByteVector &operator<<(T const &datum)
		{
			convert<T> c(datum);
			std::copy(c.raw(), c.raw() + sizeof(T), std::back_inserter(d_data));
			return *this;
		}

		template <typename T, int rank>
		ByteVector &operator<<(mVector<T, rank> const &v)
		{
			std::for_each(v.begin(), v.end(), [this] (T a)
			{
				(*this) << a;
			});
			return *this;
		}

		ByteVector &operator<<(char const *c)
		{
			(*this) << std::string(c);
			return *this;
		}

		ByteVector &operator<<(std::string const &str)
		{
			std::copy(str.begin(), str.end(), std::back_inserter(d_data));
			return *this;
		}

		template <typename T>
		ByteVector &operator>>(T &datum)
		{
			char c[sizeof(T)];

			if (unsigned(d_data.end() - d_ptr) >= sizeof(T))
			{
				std::copy(d_ptr, d_ptr + sizeof(T), c);
				datum = convert<T>(c).datum();
				d_ptr += sizeof(T);
			}

			return *this;
		}

		template <typename T, int rank>
		ByteVector &operator>>(mVector<T, rank> &v)
		{
			std::generate(v.begin(), v.end(), [this] ()
			{
				T a = 0;
				(*this) >> a;
				return a;
			});

			return *this;
		}

		ByteVector &operator>>(std::string &str)
		{
			std::string s(d_data.begin(), d_data.end());
			str = s;
			return *this;
		}

		void clear()
		{
			d_data.clear();
		}

		void to_file(std::ostream &fo)
		{
			unsigned byte_size = d_data.size();
			fo.write(reinterpret_cast<char *>(&byte_size), sizeof(unsigned));

			std::for_each(d_data.begin(), d_data.end(), [&fo] (char const &c)
			{
				fo.write(&c, 1);
			});

			fo.write(reinterpret_cast<char *>(&byte_size), sizeof(unsigned));
		}

		void from_file(std::istream &fi)
		{
			clear();

			unsigned bs1, bs2;	
			fi.read(reinterpret_cast<char *>(&bs1), sizeof(unsigned));

			std::generate_n(back_inserter(d_data), bs1, [&fi] ()
			{
				char c;
				fi.read(&c, 1);
				return c;
			});

			fi.read(reinterpret_cast<char *>(&bs2), sizeof(unsigned));
			if (bs1 != bs2)
				throw "The block seems to be corrupted";
			
			d_ptr = d_data.begin();
		}
};

} // namespace Conan

