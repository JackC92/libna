#ifndef NA_TYPE_TRAITS_IS_DENSE_XPR_H
#define NA_TYPE_TRAITS_IS_DENSE_XPR_H
#include <type_traits>
#include "Eigen/Core"
#include "na/type_traits/is_array_xpr.h"
#include "na/type_traits/is_matrix_xpr.h"

namespace na
{
	template <typename T>
	constexpr bool is_dense_xpr_v = na::is_array_xpr_v<T> || na::is_matrix_xpr_v<T>;

	template <typename T>
	struct is_dense_xpr : std::bool_constant<is_dense_xpr_v<T>> {};
}

#endif // !NA_TYPE_TRAITS_IS_DENSE_XPR_H
