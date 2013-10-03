#pragma once

#include <iostream>

namespace Misc
{
	inline std::string center(std::string const &s, int n)
	{
		int l = s.length();
		std::ostringstream o;
		o 	<< std::string((n - l) / 2, ' ') << s 
			<< std::string((n - (n - l)/2 - l), ' ');
		return o.str();
	}

	class ProgressBar
	{
		size_t N, i, j;
		std::string txt;

		public:
			ProgressBar(int N_, std::string const &txt_ = "crunching ..."):
				N(N_), i(0), j(0), txt(center(txt_, 50))
			{
				draw();
			}

			void tic()
			{
				++i;
				if (i > N) i = N;
				else while (i * 50 / N > j)
				{
					++j;
					draw();
				}
			}

			void finish()
			{
				j = 50;
				draw();
				std::cerr << std::endl;
			}

			void draw() const
			{
				std::cerr << "\r\033[m(\033[44;33;1m"
					<< txt.substr(0, j) << "\033[m"
					<< txt.substr(j, 50) << ")";
			}
	};

	template <typename Iter>
	class Progress: public Iter
	{
		ProgressBar PB;

		public:
			Progress(Iter const &i, int N, std::string const &txt): 
				Iter(i), PB(N, txt) {}

			Progress &operator++() { Iter::operator++(); PB.tic(); return *this; }

			~Progress() { PB.finish(); }
	};
}

