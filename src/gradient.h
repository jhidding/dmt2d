#pragma once

#include "conan/mvector.h"
#include "conan/cube.h"
#include "conan/boxconfig.h"
#include "cvector/cvector.h"

template <int rank>
class Gradient: public Conan::Cube<Conan::mVector<double, rank>> 
{
	public: 
		Gradient(Conan::Cube<double> const &f)
		{
			cVector dx[rank];
			for (unsigned k = 0; k < rank; ++k)
			{
				dx[k] = cVector(rank, 1, 1 << k);
			}

			size_t size = Conan::BoxConfig::size();
			unsigned bits = Conan::BoxConfig::bits();
			for (cVector x(rank, bits, 0); x < size; ++x)
			{
				for (unsigned k = 0; k < rank; ++k)
				{
					(*this)[x][k] = (f[x + dx[k]] - f[x - dx[k]]) / 2 / Conan::BoxConfig::scale();
				}
			}
		}
};

