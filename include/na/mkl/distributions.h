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

		enum Distribution
		{
			DIST_UNIFORM,
			DIST_GAUSSIAN,
			DIST_GAUSSIANMV,
			DIST_EXPONENTIAL,
			DIST_LAPLACE,
			DIST_WEIBULL,
			DIST_CAUCHY,
			DIST_RAYLEIGH,
			DIST_LOGNORMAL,
			DIST_GUMBEL,
			DIST_GAMMA,
			DIST_BETA,
			DIST_CHISQUARE
		};

		template <typename Scalar, int Distribution>
		struct distribution_generator {};

		NA_MAKE_MKL_GENERATOR(int,    DIST_UNIFORM, i, Uniform, const int, const int);
		NA_MAKE_MKL_GENERATOR(float,  DIST_UNIFORM, s, Uniform, const float, const float);
		NA_MAKE_MKL_GENERATOR(double, DIST_UNIFORM, d, Uniform, const double, const double);
		
#undef NA_MAKE_MKL_GENERATOR
	}
}

#endif // !defined(NA_MKL_DISTRIBUTIONS_H) && defined(NA_USE_MKL)
