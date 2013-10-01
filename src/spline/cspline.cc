#include "interpol.h"

double Spline::Interpol<1>::A[16] = {
	 1,  0,  0,  0,
	 0,  0,  1,  0,
	-3,  3, -2, -1,
	 2, -2,  1,  1
};

