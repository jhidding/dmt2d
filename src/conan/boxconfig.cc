#include "boxconfig.h"

// implementation ========================================================

using namespace Conan;

BoxConfig *BoxConfig::singleton = NULL;

BoxConfig::BoxConfig(unsigned rank, unsigned N, float L):
	d_N(N), d_bits(log_binary(N)), d_rank(rank),
	d_size(pow_binary(1, rank * d_bits)),
	d_L(L), d_M(rank, d_bits), d_scale(L/N), 
	d_scale2(d_scale * d_scale)
{}

BoxConfig const &BoxConfig::initialize(unsigned rank_, unsigned N_, float L_)
{
	if (singleton == NULL)
	{
		singleton = new BoxConfig(rank_, N_, L_);
		return *singleton;
	}
	else
	{
		throw;
	}
}

BoxConfig const &BoxConfig::instance()
{
	return *singleton;
}
