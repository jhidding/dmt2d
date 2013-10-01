#pragma once

inline int ipow(int base, int exponent)
{
	if (exponent == 0) return 1;
	int a = base;
	while (--exponent > 0)
	{
		a *= base;
	}
	return a;
}

