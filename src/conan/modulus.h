#pragma once

namespace Conan
{

/*!
 * calculates a mod b, by the mathematical definition:
 * the C operator % does not give the right answer if a < 0
 */
inline int modulus(int a, int b)
{
	if (a < 0) return b + (a %  b);
	else return (a % b);
}

/*!
 * calculates the a mod b, where a and b are double floats.
 */
inline double modulus(double a, double b)
{
	return a - floor(a / b) * b;
}

} // namespace Conan

