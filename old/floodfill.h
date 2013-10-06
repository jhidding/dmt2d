#pragma once

#include <queue>
#include <vector>
#include "cvector/cvector.h"

namespace Conan {

/*!
 * \brief generic class to perform flood fill
 *
 * The functionality of the flood fill should be implemented in the predicate
 * that is given with the unary_fill() call. This predicate object should have
 * overloaded the call operator to return a bool.
 * When this functor return true, the FloodFill object will try its neigbours
 * as well; if returned false, the FloodFill stops there, returning to some previous
 * point (using the std::queue structure).
 * The call to the predicate might trigger all kinds of side effects, this can be
 * a simple classic flood fill, but also more complicated behaviour.
 */
template <int rank>
class FloodFill
{
	std::vector<cVector> nb;

	public:
		FloodFill(unsigned connectivity = 4)
		{
			unsigned bits = BoxConfig::bits();

			switch (connectivity)
			{
				case 4: make_neighbours_4(nb, rank, bits);
					break;
				case 8: make_neighbours_8(nb, rank, bits);
					break;
				default: throw;
			}
		}

		template <typename UnaryPredicate>
		bool unary_fill(UnaryPredicate pred, cVector start)
		{
			std::queue<cVector> 	Q;

			if (pred(start))
			{
				Q.push(start);
			}
			else
			{
				return false;
			}

			while (!Q.empty())
			{
				cVector c = Q.front();
				
				for (unsigned i = 1; i < nb.size(); ++i)
				{
					cVector p = c + nb[i];

					if (pred(p))
					{
						Q.push(p);
					}
				}

				Q.pop();
			}

			return true;
		}
};

} // namespace Conan

