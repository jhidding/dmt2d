#include "interpol.h"
#include <iostream>

using namespace std;
using namespace Spline;
using namespace Conan;

Interpol<3>::Cache::Cache(Cache const &o):
	f(o.f)
{
	cache_size = 4096;
	cache_data = new double[cache_size * 4096];
	cache_ptr = cache_data;
}

Interpol<3>::Cache::Cache(double const *f_): 
	f(f_) 
{
	cache_size = 4096;
	cache_data = new double[cache_size * 4096];
	cache_ptr = cache_data;
}

double *Interpol<3>::Cache::operator[](cVector const &idx)
{
	if (C.count(idx) == 1)
	{
		return C[idx];
	}
	else
	{
		if (Q.size() >= cache_size)
		{
			cache_ptr = pop();
		}

		if (static_cast<unsigned>(cache_ptr - cache_data) > ((cache_size - 1) * 4096))
		{
			cerr << Q.size() << " " << cache_ptr - cache_data << " " << cache_size << endl;
			throw "weird memory overflow";
		}

		push(idx);
		cache_ptr += 4096;

		return C[idx];
	}
}

void Interpol<3>::Cache::push(cVector const &S)
{
	cVector dx(3, 1, 1), dy(3, 1, 2), dz(3, 1, 4);
	double b[64];

	for (cVector i(3, 1, 0); i < 8; ++i)
	{
		cVector X = S + i;

		// f
		b[i]		= f[X];

		// f_x
		b[i + 8] 	= (f[X + dx] - f[X - dx]) / 2;

		// f_y
		b[i + 16] 	= (f[X + dy] - f[X - dy]) / 2;

		// f_z
		b[i + 24] 	= (f[X + dz] - f[X - dz]) / 2;

		// f_xy
		b[i + 32]	= (f[X + dx + dy] - f[X + dx - dy]
				-  f[X - dx + dy] + f[X - dx - dy]) / 4;

		// f_yz
		b[i + 40]	= (f[X + dy + dz] - f[X + dy - dz]
				-  f[X - dy + dz] + f[X - dy - dz]) / 4;

		// f_zx
		b[i + 48]	= (f[X + dz + dx] - f[X + dz - dx]
				-  f[X - dz + dx] + f[X - dz - dx]) / 4;

		// f_xyz
		b[i + 56]	= (f[X + dx + dy + dz] - f[X + dx + dy - dz]
				-  f[X + dx - dy + dz] + f[X + dx - dy - dz]
				-  f[X - dx + dy + dz] + f[X - dx + dy - dz]
				+  f[X - dx - dy + dz] - f[X - dx - dy - dz]) / 8;
	}

	double *alpha = cache_ptr;
	std::fill_n(alpha, 64, 0.0);
	for (cVector i(2, 6, 0); i < 4096; ++i)
		alpha[i[1]] += b[i[0]] * A[i];

	C[S] = alpha;
	Q.push(S);
}

double *Interpol<3>::Cache::pop()
{
	double *a = C[Q.front()];
	C.erase(Q.front());
	Q.pop();
	return a;
}

Interpol<3>::Cache::~Cache()
{
	delete[] cache_data;
}

Interpol<3>::Interpol(double const *f_): bits(Conan::BoxConfig::bits()), cache(f_) {}

double Interpol<3>::operator()(mVector<double, rank> const &x)
{
	return f(x);
}

double Interpol<3>::f(mVector<double, rank> const &x)
{
	mVector<double, rank> P;
	std::transform(x.begin(), x.end(), P.begin(), fraction);

	mVector<int, rank> start;
	std::transform(x.begin(), x.end(), start.begin(), integer);
	cVector S(rank, bits, BoxConfig::M(start));

	double *alpha = cache[S];

	// P(x, y, z) = Sum_i Sum_j Sum_k [alpha_ijk * x^i y^j z^k]
	double v = 0;
	for (cVector i(rank, 2, 0); i < (1 << (2*rank)); ++i) // convenient way to do triple loop
		v += alpha[i] * pow(P[0], i[0]) * pow(P[1], i[1]) * pow(P[2], i[2]);

	return v;
}

double Interpol<3>::df(unsigned k, mVector<double, rank> const &x)
{
	mVector<double, rank> P;
	std::transform(x.begin(), x.end(), P.begin(), fraction);

	mVector<int, rank> start;
	std::transform(x.begin(), x.end(), start.begin(), integer);
	cVector S(rank, bits, BoxConfig::M(start));

	double *alpha = cache[S];

	double v = 0;
	unsigned i = (k + 1) % 3;
	unsigned j = (k + 2) % 3;

	for (cVector z(rank, 2, 0); z < (1 << (2 * rank)); ++z) // convenient way to do triple loop
		if (z[k] != 0)
			v += alpha[z] * z[k] * pow(P[k], z[k] - 1) * pow(P[i], z[i]) * pow(P[j], z[j]);

	return v;
}

