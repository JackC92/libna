#ifndef NA_ADAPTIVE_SIMPSON_H
#define NA_ADAPTIVE_SIMPSON_H
#include <cmath>
#include <complex>
#include <functional>
#include "na/type_traits/floating_point_scalar.h"
#include "na/type_traits/is_floating_point_or_complex.h"

namespace na
{
	namespace internal
	{
		template <typename Scalar>
		Scalar adaptive_simpson_aux(
			const std::function<Scalar(double)> f,
			const double a,
			const double b,
			const double tol,
			const Scalar whole,
			const Scalar fa,
			const Scalar fb,
			const Scalar fm,
			const int depth)
		{
			typedef typename na::floating_point_scalar<Scalar>::value_type Real;
			double m = 0.5 * (a + b), h = 0.5 * (b - a);
			double lm = 0.5 * (a + m), rm = 0.5 * (m + b);
			Scalar flm = f(lm), frm = f(rm);
			Scalar left = (h / 6.0) * (fa + 4.0 * flm + fm);
			Scalar right = (h / 6.0) * (fm + 4.0 * frm + fb);
			Scalar delta = (left + right) - whole;
			if ((depth <= 0) || (std::abs(delta) <= 15.0 * tol))
			{
				return left + right + (delta / static_cast<Real>(15.0));
			}
			else
			{
				return adaptive_simpson_aux(f, a, m, 0.5 * tol, left, fa, fm, flm, depth - 1)
					+ adaptive_simpson_aux(f, m, b, 0.5 * tol, right, fm, fb, frm, depth - 1);
			}
		}
	}
	
	template <typename Scalar>
	Scalar adaptive_simpson(
		const std::function<Scalar(double)> f,
		const double a,
		const double b,
		const double tol,
		const int max_depth)
	{
		static_assert(na::is_floating_point_or_complex_v<Scalar>, "Scalar must be a floating-point type or std::complex");
		double h = b - a;
		if (h == 0.0)
		{
			return Scalar();
		}
		Scalar fa = f(a), fb = f(b), fm = f(0.5 * (a + b));
		Scalar whole = (h / 6.0) * (fa + 4.0 * fm + fb);
		return na::internal::adaptive_simpson_aux(f, a, b, tol, whole, fa, fb, fm, max_depth);
	}
}

#endif // !NA_ADAPTIVE_SIMPSON_H
