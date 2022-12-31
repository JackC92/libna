#if !defined(NA_MKL_DISTRIBUTIONS_H) && defined(NA_USE_MKL)
#define NA_MKL_DISTRIBUTIONS_H
#include "mkl_vsl.h"
#include "na/macros.h"

namespace na
{
	namespace mkl
	{

#define NA_MAKE_MKL_GENERATOR(Type, Distribution, Symbol, Suffix, ...) \
EXPSPEC struct distribution_generator<Type, Distribution> \
{ \
	typedef int (*GeneratorPtr)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, Type[], __VA_ARGS__); \
	static inline GeneratorPtr run = &v ## Symbol ## Rng ## Suffix; \
};

		enum class Distribution
		{
			UNIFORM,
			GAUSSIAN,
			GAUSSIANMV,
			EXPONENTIAL,
			LAPLACE,
			WEIBULL,
			CAUCHY,
			RAYLEIGH,
			LOGNORMAL,
			GUMBEL,
			GAMMA,
			BETA,
			CHISQUARE
		};

		template <typename Scalar, Distribution dist>
		struct distribution_generator {};

		NA_MAKE_MKL_GENERATOR(int,    Distribution::UNIFORM, i, Uniform, const int, const int);
		NA_MAKE_MKL_GENERATOR(float,  Distribution::UNIFORM, s, Uniform, const float, const float);
		NA_MAKE_MKL_GENERATOR(double, Distribution::UNIFORM, d, Uniform, const double, const double);
		
#undef NA_MAKE_MKL_GENERATOR
	}
}

#endif // !defined(NA_MKL_DISTRIBUTIONS_H) && defined(NA_USE_MKL)
