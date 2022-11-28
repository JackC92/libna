#ifndef NA_SPECIAL_BESSEL_H
#define NA_SPECIAL_BESSEL_H

// Following the notation of Mathematica, we define
//  - BesselJ[n, x]: the nth-order Bessel function of the first kind J_n(x), where n is a non-negative integer and x is a non-negative real number.
//  - BesselJZero[n, k]: the kth zero of the Bessel function of the first kind J_n(x), where n is a non-negative integer and k is a positive integer.
namespace na
{
	namespace bessel
	{
		// Compute the value BesselJZero[0, k]
		//
		// Input parameters:
		//   k: the index of the zero of BesselJ[0, x]
		//
		// Output parameters:
		//   val: the value of the zero of BesselJ[0, x]
		double besselj0_zero(const int k);

		// Compute the square of BesselJ[1, BesselJZero[0, k]]
		//
		// Input parameters:
		//   k: the index of the zero of BesselJ[0, x]
		//
		// Output parameters:
		//   val: the square of BesselJ[1, BesselJZero[0, k]]
		double besselj1_squared(const int k);
	}
}

#endif // !NA_SPECIAL_BESSEL_H
