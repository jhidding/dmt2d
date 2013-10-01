#pragma once

#include <complex>

#ifndef __complex64__
#define __complex64__
typedef std::complex<double> complex64;
#endif

extern std::complex<double> math_i;
inline double get_real(std::complex<double> &z) { return z.real(); }
inline double get_imag(std::complex<double> &z) { return z.imag(); }

