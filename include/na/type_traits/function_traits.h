#ifndef NA_TYPE_TRAITS_FUNCTION_TRAITS_H
#define NA_TYPE_TRAITS_FUNCTION_TRAITS_H
#include <tuple>

namespace na
{
	template <typename T>
	struct function_traits;

	template <typename T>
	struct function_traits<const T> : function_traits<T> {};

	template <typename R, typename... Args>
	struct function_traits<std::function<R(Args...)>>
	{
		static constexpr std::size_t arity = sizeof...(Args);

		using return_type = R;
		using arg_types = std::tuple<Args...>;

		template <std::size_t N>
		struct arg
		{
			static_assert(N < arity, "arg: invalid parameter index.");
			using type = typename std::tuple_element<N, arg_types>::type;
		};
	};
}

#endif // !NA_TYPE_TRAITS_FUNCTION_TRAITS_H
