#pragma once

#include "bit_tricks.h"

#include <map>
// #include <iostream>

/*
 * cVector_info stores information on a certain class of
 * cVectors. This information is too much to store in every
 * instance of a cVector, but changes if the bit-paramater
 * is different.
 *
 * cVector_info_map maps for each rank the cVector_info
 * to some bit-parameter. When a cVector is initialized
 * it looks up which information to use and saves a pointer
 * to that information.
 */

class cVector
{
	class Info
	{
		public:
			BitTricks::ReverseBitorder rev;
			BitTricks::MortonExpansion mor;
			BitTricks::InverseLog ilg;

			unsigned	*bits;
			std::size_t		*mask;
			std::size_t		size;

			Info(unsigned rank, unsigned b);
			~Info();

		private:
			Info(Info const &); // NI
			Info &operator=(Info const &); // NI
	};

	class InfoMap: public std::map<unsigned, cVector::Info *>
	{
		public:
			InfoMap():
				std::map<unsigned, cVector::Info*>() {}
			Info const *operator[](unsigned idx);
			~InfoMap();
	};

	static InfoMap	info_map;
	Info const 	*info;

	std::size_t		X;

	public:
		unsigned	bits;
		unsigned 	rank;

		cVector() {}
		cVector(unsigned r, unsigned b_, std::size_t X_ = 0);

		// get one element of the vector
		inline unsigned operator[](unsigned idx) const;

		// increment and decrement operators
		// if there is a cary it will add up to a
		// higher dimension, so you can traverse an
		// entire array in this fashion.
		cVector &operator--() { --X; return *this; }
		cVector &operator++() { ++X; return *this; }
		cVector operator--(int) { cVector a(*this); --X; return a; }
		cVector operator++(int) { cVector a(*this); ++X; return a; }

		// add or substract two vectors that can have
		// different bits but must have same rank
		cVector &operator+=(cVector const &o);
		cVector &operator-=(cVector const &o);

		// copying versions of addition and substraction
		cVector operator+(cVector const &o) const
		{ return cVector(*this) += o; }

		cVector operator-(cVector const &o) const
		{ return cVector(*this) -= o; }

		// these comparisons are only here to suppress compiler
		// warnings on dubious conversion. Just remember that
		// a cVector is a std::size_t with some strings attached.
		bool operator!=(unsigned s) const { return X != s; }
		bool operator==(unsigned s) const { return X == s; }

		bool operator!=(cVector const &s) const { return X != s.X; }
		bool operator==(cVector const &s) const { return X == s.X; }

		bool operator<(cVector const &s) const { return X < s.X; }
		bool operator>(cVector const &s) const { return X > s.X; }
		bool operator<=(cVector const &s) const { return X <= s.X; }
		bool operator>=(cVector const &s) const { return X >= s.X; }

		operator std::size_t() const { return X; }

		std::size_t size() const;

		cVector reverse() const;
		cVector reverse_morton() const;
		cVector morton() const;
		cVector invlog() const;
};

inline unsigned cVector::operator[](unsigned idx) const
{ return (X >> (bits * idx)) & ((1 << bits) - 1); }
//		& info->mask[idx]) >> info->bits[idx]; }

extern void make_neighbours_4(std::vector<cVector> &N, unsigned rank, unsigned bits);
extern void make_neighbours_8(std::vector<cVector> &N, unsigned rank, unsigned bits);

