#include "na/special/gamma.h"
#include <cassert>
#include <cmath>

namespace na
{
	namespace gamma
	{
		double binomial_coefficient(
			const int n,
			const int k)
		{
			assert((n >= 0) && "binomial_coefficient: n must be non-negative");
			assert((k >= 0) && "binomial_coefficient: k must be non-negative");

			return std::floor(0.5 + std::exp(std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0)));
		}
	}
}
