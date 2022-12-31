#ifndef NA_TYPE_TRAITS_IS_POINTER_LIKE_H
#define NA_TYPE_TRAITS_IS_POINTER_LIKE_H
#include <type_traits>
#include "na/type_traits/is_detected.h"

namespace na
{
	template <typename T>
	using deferenced_type = decltype(*std::declval<T>());

	template <typename T>
	using arrow_type = decltype(std::declval<T>().operator->());
	
	template <typename T>
	using pointed_type = detected_or_t<detected_t<deferenced_type, T>, arrow_type, T>;

	template <typename T>
	constexpr bool is_pointer_like_v = std::disjunction_v<
		na::is_detected<deferenced_type, T>,
		na::is_detected<arrow_type, T>>;

	template <typename T>
	struct is_pointer_like : std::bool_constant<is_pointer_like_v<T>> {};
}

#endif // !NA_TYPE_TRAITS_IS_POINTER_LIKE_H
