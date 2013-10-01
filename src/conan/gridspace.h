#pragma once

//#include "../conan.h"
#include "ndspace.h"
#include "map.h"
#include "mvector.h"

//namespace Spaces {
namespace Conan {

template <typename T, int rank>
class GridSpace: public Map<NdSpace<rank>, mVector<T, rank>>
{
	NdSpace<rank> S;
	public:
		GridSpace(unsigned N, float L)
			: Map<NdSpace<rank>, mVector<T, rank>>(
				[N, L] (mVector<T, rank> const &k) { return k * (L / N); },
				S),
			S(mVector<int, rank>(N))
		{}
};
}
//} // namespace Spaces

