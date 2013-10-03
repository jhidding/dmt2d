#pragma once
#include "system.hh"

namespace Misc
{
	template <typename Q, unsigned R>
	class Interpolate
	{
		using dVector = System::mVector<double, R>;
		using iVector = System::mVector<int, R>;

		System::ptr<System::BoxConfig<R>> 	box;
		Q									f;

		public:
			Interpolate(System::ptr<System::BoxConfig<R>> box_, Q f_):
				box(box_), f(f_)
			{}

			typename Q::value_type operator()(dVector const &x) const
			{
				iVector origin;
				dVector A[2];
				double r = box->scale();

				transform(x, origin, [r] (double q)
					{ return (int)(q); });
				transform(x, origin, A[1], [] (double q, int o) 
					{ return q - o; });
				transform(A[1], A[0], [] (double g)
					{ return 1 - g; });

				auto s = box->m2c(origin);
				typename Q::value_type v(0);

				for (auto dx : box->sq())
				{
					double z = 1;
					for (unsigned k = 0; k < R; ++k) 
						z *= A[dx[k]][k];

					v += f[s + dx] * z;
				}

				return v;
			}
	};
}

// vim:sw=4:ts=4:tw=72
