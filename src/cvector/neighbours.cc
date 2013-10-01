#include "cvector.h"
#include <vector>
#include <algorithm>

#include "ipow.h"

/*
 * connectivity of cvectors
 *
 * allways the first cvector in &N will be the centre cell.
 *
 * 4: get the perpendicular ajacent cells. In 2-D there are 4, in 3-D 6.
 * 8: is more complicated, we need 8 neighbours (including the diagonals)
 * in 2-D and 26 in 3-D. Each cell starts in the middle. Then for each direction
 * we take the top third that gets +1, the middle third gets +0 and the bottom
 * third gets -1.
 */

void make_neighbours_4(std::vector<cVector> &N, unsigned rank, unsigned bits)
{
	N.push_back(cVector(rank,bits,0));
	for (unsigned i = 0; i < rank; ++i)
	{
		N.push_back(cVector(rank, bits, std::size_t(1) << (i * bits)));
		N.push_back(cVector(rank, bits, ((std::size_t(1) << bits) - 1) << (i * bits)));
	}
}

// comment shows 3-D case
void make_neighbours_8(std::vector<cVector> &N, unsigned rank, unsigned bits)
{
	std::size_t v[3];
	// +0
	v[0] = 0;
	// +1
	v[1] = std::size_t(1);
	// -1
	v[2] = (std::size_t(1) << bits) - 1;

	// 27 cells
	N.resize(ipow(3, rank));
	// fill with 0
	std::fill(N.begin(), N.end(), cVector(rank, bits, 0));

	// loop over directions ...
	for (unsigned i = 0; i < rank; ++i)
	{
		unsigned j = 0;
		// loop l = 0 ... 1, 3, 9
		for (int l = 0; l < ipow(3, i); ++l)
		{
			// loop k = 0 ... 3
			for (unsigned k = 0; k < 3; ++k)
			{
				// loop m = 0 ... 9, 3, 1
				for (int m = 0; m < ipow(3, rank - i - 1); ++m)
				{
					N[j] += cVector(rank, bits, v[k] << (i * bits));
					++j;
				}
			}
		}
	}
}

