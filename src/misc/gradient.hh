#pragma once
#include "../base/mvector.hh"
#include "../base/boxconfig.hh"

namespace Misc
{
	template <unsigned R>
	class Gradient
	{
		typedef System::mVector<double, R> Vector;

		System::Array<double> data;
		System::ptr<System::cVector<R>> b;

		public:
			typedef Vector value_type;

			Gradient(System::ptr<System::BoxConfig<R>> box_, System::Array<double> data_):
				data(data_), b(new System::cVector<R>(box_->bits())) {}

			inline double fdi(size_t i, unsigned k) const
			{
				return data[b->sub(i, b->dx2[k])] / 12. - 2. * data[b->sub(i, b->dx[k])] / 3.
					+ 2. * data[b->add(i, b->dx[k])] / 3. - data[b->add(i, b->dx2[k])] / 12.;
			}

			inline Vector grad(size_t i) const
			{
				Vector v; 
				for (unsigned k = 0; k < R; ++k)
					v[k] = fdi(i, k);
				return v;
			}

			bool root_potentially_within_cell(size_t i) const
			{
				// true if within the cell all components of
				// the gradient change sign.
				std::vector<bool> test(R, false);
				int cnt = 0;

				Vector v1 = grad(i);
				for (auto &dx : b->sq)
				{
					Vector v2 = grad(b->add(i, dx));
					for (unsigned k = 0; k < R; ++k)
					{
						if (v1[k] * v2[k] < 0. and not test[k])
						{
							test[k] = true;
							++cnt;
						}

						if (cnt == R) return true;
					}
				}

				return false;
			}

			Vector find_root(size_t i) const
			{
				return Vector();
			}

			Vector operator[](size_t i) const
			{
				return grad(i);
			}
	};
}

