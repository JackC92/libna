#ifndef NA_TYPE_TRAITS_IS_MATRIX_XPR_H
#define NA_TYPE_TRAITS_IS_MATRIX_XPR_H
#include <type_traits>
#include "Eigen/Core"

namespace na
{
	template <typename T>
	constexpr bool is_matrix_xpr_v = false;

	template <typename T>
	constexpr bool is_matrix_xpr_v<Eigen::MatrixBase<T>> = true;

	template <typename T>
	struct is_matrix_xpr : std::bool_constant<is_matrix_xpr_v<T>> {};
}

#endif // !NA_TYPE_TRAITS_IS_MATRIX_XPR_H
