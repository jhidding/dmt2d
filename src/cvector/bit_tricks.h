#pragma once

#include <vector>

namespace BitTricks
{

template <typename T>
class Counter
{
	T value;
	public:
		Counter(T init): value(init) {}
		T operator()() { return value++; }
};

/*
 * ReverseBitorder creates a lookup table for the reversion
 * of the bits in a number. The number of bits that are reversed
 * can be arbitrary but should not be too big or bigger than 
 * log2(max_int) - 1. It is intended to be used as a tool to perform
 * an iterative binary grid walk going 0 1/2 1/4 3/4 1/8 3/8 5/8 etc.
 */

class ReverseBitorder: public std::vector<std::size_t>
{
	unsigned bits;
	public:
		ReverseBitorder(unsigned b);
		std::size_t operator()(std::size_t i) const;
};

/*
 * Morton numbers are a fancy way to index on a grid.
 * This class builds a lookup table to get a
 * "Morton expansion" of a number. That is:
 * 	000.000.000.abc -> 000a.000b.000c
 *
 * So we need two parameters: rank and bits
 */

class MortonExpansion: public std::vector<std::size_t>
{
	unsigned rank, bits;

	public:
		MortonExpansion(unsigned r, unsigned b);
		std::size_t operator()(std::size_t i);
};

class InverseLog: public std::vector<std::size_t>
{
	unsigned bits;
	
	public:
		InverseLog(unsigned b);
		std::size_t operator()(std::size_t i) const;
};

} // BitTricks

