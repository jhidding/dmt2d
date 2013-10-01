#pragma once

union cVector_type_union
{
	struct
	{
		unsigned char rank, bits;
	} p;

	unsigned value;
};

inline unsigned cVector_type(unsigned char rank, unsigned char bits)
{
	cVector_type_union U;
	U.p.rank = rank;
	U.p.bits = bits;

	return U.value;
}

