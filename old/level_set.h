#pragma once
#include "conan/mvector.h"
#include "conan/boxconfig.h"
#include <iostream>
#include <vector>

class Contour_2: public std::vector<Conan::mVector<double, 2>>
{
	typedef Conan::mVector<double, 2> Point;
	
	public:
		void output_for_gnuplot(std::ostream &out) const
		{
			for (unsigned a = 0; a <= size(); ++a)
			{
				unsigned i = Conan::modulus(int(a), int(size())),
					 j = Conan::modulus(int(a) - 1, int(size()));

				if (Conan::norm((*this)[i] - (*this)[j]) > Conan::BoxConfig::L() / 2)
					out << "\n\n";

				out << (*this)[i] << "\n";
			}
		}

		template <typename F>
		Contour_2 map(F f) const
		{
			Contour_2 result;
			std::transform(begin(), end(), std::back_inserter(result), f);
			return result;
		}
};

class Level_set_2: public std::vector<Contour_2>
{
	typedef Conan::mVector<double, 2> Point;

	public:
		Level_set_2(double const *f);
		void output_for_gnuplot(std::ostream &out) const;
};

class Double_level_set_2: public std::vector<Conan::mVector<double, 2>>
{
	typedef Conan::mVector<double, 2> Point;

	public:	
		Double_level_set_2() {}
		Double_level_set_2(double const *f, double const *g);

		template <typename F>
		Double_level_set_2 map(F f) const
		{
			Double_level_set_2 result;
			std::transform(begin(), end(), std::back_inserter(result), f);
			return result;
		}

		void output_for_gnuplot(double const *ev, std::ostream &out) const;
		void output_for_gnuplot(std::ostream &out) const;
};

#include "conan/interpol.h"

class Critical_line_2: 
	public std::vector<Conan::mVector<double, 2>>
{
	typedef Conan::mVector<double, 2> Point;
	typedef std::vector<Point> Line;

	enum { UNSET, IS_LOOP, HAS_ENDS } line_type;
	std::pair<Point, Point> ends;
	bool start_set, end_set;

	public:
		Critical_line_2(): line_type(UNSET), 
			start_set(false), end_set(false) {}

		void set_to_being_loop() { line_type = IS_LOOP; }
		void set_to_have_ends() { line_type = HAS_ENDS; }
		void set_start(Point const &Q) { ends.first = Q; start_set = true; }
		void set_end(Point const &Q) { ends.second = Q; end_set = true; }
		bool is_loop() const { return line_type == IS_LOOP; }
		bool has_ends() const { return line_type == HAS_ENDS; }
		bool unset() const { return line_type == UNSET; }
		bool start_unset() const { return not start_set; }
		bool end_unset() const { return not end_set; }

		void output_for_gnuplot(double const *ev, std::ostream &out) const;
		template <typename F>
		void output_mapped(F f, double const *Eval, std::ostream &out) const
		{
			Conan::Interpol<double, 2> ev(Eval);

			size_t s = size();
			if (is_loop()) s++;

			bool pen_on_paper = true;
			for (unsigned a = 0; a < s; ++a)
			{
				unsigned i = Conan::modulus(int(a), int(size())),
					 j = Conan::modulus(int(a) - 1, int(size()));

				if (Conan::norm(f((*this)[i]) - f((*this)[j])) 
						> Conan::BoxConfig::L() / 2)
					out << "\n\n";

				double v = ev((*this)[i]);
				if (v > 0.0)
				{
					if (not pen_on_paper) out << "\n\n";
					out << f((*this)[i]) << " " << ev((*this)[i]) << "\n";
					pen_on_paper = true;
				}
				else
				{
					pen_on_paper = false;
				}
			}
		}
};

class A3_lines_2
{
	typedef Conan::mVector<double, 2> Point;
	typedef Critical_line_2 Line;
	typedef std::vector<Line> Line_set;
	
	Line_set data;

	public:
		A3_lines_2(double const *Eval,
			double const *H11, double const *H12, double const *H22,	
			std::vector<Point> &D3_points, bool, bool);

		std::vector<Point> A4_points(double const *Eval);

		void output_for_gnuplot(double const *ev, std::ostream &out) const;
		template <typename F>
		void output_mapped(F f, double const *ev, std::ostream &out) const
		{
			std::for_each(data.begin(), data.end(),
				[&] (Line const &L)
			{
				L.output_mapped(f, ev, out);
				out << "\n\n";
			});
		}
};


