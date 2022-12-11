#ifndef NA_TYPE_TRAITS_IS_ARRAY_OR_MATRIX_H
#define NA_TYPE_TRAITS_IS_ARRAY_OR_MATRIX_H
#include <type_traits>
#include "Eigen/Core"

namespace na
{
	template <typename T>
	constexpr bool is_array_or_matrix_v = false;

	template <typename T>
	constexpr bool is_array_or_matrix_v<Eigen::PlainObjectBase<T>> = true;

	template <typename T>
	struct is_array_or_matrix : std::bool_constant<is_array_or_matrix_v<T>> {};
}

#endif // !NA_TYPE_TRAITS_IS_ARRAY_OR_MATRIX_H
