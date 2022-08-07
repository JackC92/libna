#ifndef NA_TYPE_TRAITS_IS_COMPLEX_H
#define NA_TYPE_TRAITS_IS_COMPLEX_H
#include <complex>
#include <type_traits>

namespace na
{
	template <typename T>
	constexpr bool is_complex_v = std::disjunction_v<
		std::is_same<std::remove_cv_t<T>, std::complex<float>>,
		std::is_same<std::remove_cv_t<T>, std::complex<double>>,
		std::is_same<std::remove_cv_t<T>, std::complex<long double>>>;

	template <typename T>
	struct is_complex : std::bool_constant<is_complex_v<T>> {};
}

#endif // !NA_TYPE_TRAITS_IS_COMPLEX_H
