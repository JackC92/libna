#ifndef NA_TYPE_TRAITS_IS_FLOATING_POINT_OR_COMPLEX_H
#define NA_TYPE_TRAITS_IS_FLOATING_POINT_OR_COMPLEX_H
#include <complex>
#include <type_traits>
#include "na/type_traits/is_complex.h"

namespace na
{
	template <typename T>
	constexpr bool is_floating_point_or_complex_v = std::is_floating_point_v<T> || na::is_complex_v<T>;

	template <typename T>
	struct is_floating_point_or_complex : std::bool_constant<is_floating_point_or_complex_v<T>> {};
}

#endif // !NA_TYPE_TRAITS_IS_FLOATING_POINT_OR_COMPLEX_H
