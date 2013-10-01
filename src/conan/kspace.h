#pragma once

#include <cmath>
#include <algorithm>
#include "ndspace.h"
#include "map.h"
#include "mvector.h"

/*!
 * \ingroup spaces
 * Maps a CalcK<rank> to NdSpace<rank> to create the impression
 * of an array of KSpace values.
 */

namespace Conan {

template <int rank>
class KSpace: public Map<NdSpace<rank>, mVector<double, rank>>
{
	NdSpace<rank> S;
	public:
		KSpace(unsigned N, double L): 	
			Map<NdSpace<rank>, mVector<double, rank>>(
				[N, L] (mVector<int, rank> const &v) 
				{
					mVector<double, rank> w;
					std::transform(v.begin(), v.end(), w.begin(), [N, L] (int i) -> double 
					{ 
						return (i > int(N/2) ? i - int(N) : i) * (2 * M_PI / L); 
					} );

					return w;
				}, S),

				S(mVector<int, rank>(N))
		{}
};

}

