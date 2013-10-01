#include "interpol.h"
#include <iostream>

using namespace std;

using namespace Spline;
using namespace Conan;

Interpol<2>::Cache::Cache(Cache const &o):
	f(o.f)
{
	cache_size = 128;
	cache_data = new double[cache_size * 16];
	cache_ptr = cache_data;
}

Interpol<2>::Cache::Cache(double const *f_): 
	f(f_) 
{
	cache_size = 128;
	cache_data = new double[cache_size * 16];
	cache_ptr = cache_data;
}

double *Interpol<2>::Cache::operator[](cVector const &idx)
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

		if (static_cast<unsigned>(cache_ptr - cache_data) > ((cache_size - 1) * 16))
		{
			cerr << Q.size() << " " << cache_ptr - cache_data << " " << cache_size << endl;
			throw "weird memory overflow";
		}

		push(idx);
		cache_ptr += 16;

		return C[idx];
	}
}

// this is the only rank dependent method
void Interpol<2>::Cache::push(cVector const &S)
{
	cVector dx(2, 1, 1), dy(2, 1, 2);
	double b[16];

	for (cVector i(2, 1, 0); i < 4; ++i)
	{
		cVector X = S + i;
		// f(x)
		b[i] 		= f[X];
		// df/dx
		b[i + 4] 	= (f[X + dx] - f[X - dx])/2;
		// df/dy
		b[i + 8] 	= (f[X + dy] - f[X - dy])/2;
		// (d^2 f)/(dy dx)
		b[i + 12]	= ((f[(X + dy) + dx] - f[(X - dy) + dx]) - 
			(f[(X + dy) - dx] - f[(X - dy) - dx]))/4;
	}

	double *alpha = cache_ptr;
	std::fill_n(alpha, 16, 0.0);
	for (cVector i(2, 4, 0); i < 256; ++i)
		alpha[i[1]] += b[i[0]] * A[i];

	C[S] = alpha;
	Q.push(S);
}

double *Interpol<2>::Cache::pop()
{
	double *a = C[Q.front()];
	C.erase(Q.front());
	Q.pop();
	return a;
}

Interpol<2>::Cache::~Cache()
{
	delete[] cache_data;
}

Interpol<2>::Interpol(double const *f_): 
	bits(Conan::BoxConfig::bits()), cache(f_) {}

double Interpol<2>::operator()(mVector<double, 2> const &x)
{
	return f(x);
}

double Interpol<2>::f(mVector<double, 2> const &x)
{
	mVector<double, 2> P;
	std::transform(x.begin(), x.end(), P.begin(), fraction);

	mVector<int, 2> start;
	std::transform(x.begin(), x.end(), start.begin(), integer);
	cVector S(2, bits, BoxConfig::M(start));

	double *alpha = cache[S];

	// v(x, y) = Sum[alpha_ij * x**i y**j]
	double v = 0;
	for (cVector z(2, 2, 0); z < 16; ++z)
		v += alpha[z] * pow(P[0], z[0]) * pow(P[1], z[1]);

	return v;
}


double Interpol<2>::df(unsigned k, mVector<double, 2> const &x)
{
	mVector<double, 2> P;
	std::transform(x.begin(), x.end(), P.begin(), fraction);

	mVector<int, 2> start;
	std::transform(x.begin(), x.end(), start.begin(), integer);
	cVector S(2, bits, BoxConfig::M(start));

	double *alpha = cache[S];

	// v(x, y) = Sum[alpha_ij * x**i y**j]
	double v = 0;
	unsigned i = (k + 1) % 2;
	for (cVector z(2, 2, 0); z < 16; ++z)
		if (z[k] != 0)
			v += alpha[z] * z[k] *  pow(P[k], z[k] - 1) * pow(P[i], z[i]);

	return v;
}

