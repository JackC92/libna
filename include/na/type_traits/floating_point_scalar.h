#ifndef NA_TYPE_TRAITS_FLOATING_POINT_SCALAR_H
#define NA_TYPE_TRAITS_FLOATING_POINT_SCALAR_H
#include <complex>
#include <type_traits>
#include "na/type_traits/is_float_or_complex.h"

namespace na
{
	template <typename T>
	struct floating_point_scalar
	{
		static_assert(na::is_float_or_complex_v<T>, "floating_point_scalar: T must be a floating-point type or std::complex");
		using value_type = T;
	};

	template <typename T>
	struct floating_point_scalar<std::complex<T>>
	{
		static_assert(std::is_floating_point_v<T>, "floating_point_scalar: T must be a floating-point type or std::complex");
		using value_type = T;
	};
}

#endif // !NA_TYPE_TRAITS_FLOATING_POINT_SCALAR_H
