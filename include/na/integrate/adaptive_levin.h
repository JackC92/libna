#ifndef NA_INTEGRATE_ADAPTIVE_LEVIN_H
#define NA_INTEGRATE_ADAPTIVE_LEVIN_H
#include <complex>
#include <functional>
#include "Eigen/Core"

namespace na
{
	// Compute the integral of f(x) * exp(g(x)) over the interval [a, b] \in \R with
		// the user-provided tolerance tol using the Chebyshev spectral method of order k.
	double adaptive_levin(
		const std::function<double(double)> f,
		const std::function<double(double)> g,
		const double a,
		const double b,
		const double tol,
		const Eigen::Index k);

	// Compute the integral of f(x) * exp(i * g(x)) over the interval [a, b] \in \R with
	// the user-provided tolerance tol using the Chebyshev spectral method of order k.
	std::complex<double> cadaptive_levin(
		const std::function<double(double)> f,
		const std::function<double(double)> g,
		const double a,
		const double b,
		const double tol,
		const Eigen::Index k);
}

#endif // !NA_INTEGRATE_ADAPTIVE_LEVIN_H
