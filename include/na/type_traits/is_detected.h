#ifndef NA_TYPE_TRAITS_IS_DETECTED_H
#define NA_TYPE_TRAITS_IS_DETECTED_H
#include <type_traits>
#include "na/type_traits/nonesuch.h"

namespace na
{
	namespace detail
	{
		template <typename Default, typename AlwaysVoid, template <typename...> class Op, typename... Args>
		struct detector
		{
			using value_t = std::false_type;
			using type = Default;
		};

		template <typename Default, template <typename...> class Op, typename... Args>
		struct detector<Default, std::void_t<Op<Args...>>, Op, Args...>
		{
			using value_t = std::true_type;
			using type = Op<Args...>;
		};
	}

	template <template <typename...> class Op, typename... Args>
	using is_detected = typename na::detail::detector<nonesuch, void, Op, Args...>::value_t;

	template <template <typename...> class Op, typename... Args>
	constexpr bool is_detected_v = is_detected<Op, Args...>::value;
	
	template <template <typename...> class Op, typename... Args>
	using detected_t = typename na::detail::detector<nonesuch, void, Op, Args...>::type;

	template <typename Default, template <typename...> class Op, typename... Args>
	using detected_or = na::detail::detector<Default, void, Op, Args...>;

	template <typename Default, template <typename...> class Op, typename... Args>
	using detected_or_t = typename detected_or<Default, Op, Args...>::type;
	
	template <typename Expected, template <typename...> class Op, typename... Args>
	using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;

	template <typename Expected, template <typename...> class Op, typename... Args>
	constexpr bool is_detected_exact_v = is_detected_exact<Expected, Op, Args...>::value;

	template <typename To, template <typename...> class Op, typename... Args>
	using is_detected_convertible = std::is_convertible<detected_t<Op, Args...>, To>;
	
	template <typename To, template <typename...> class Op, typename... Args>
	constexpr bool is_detected_convertible_v = is_detected_convertible<To, Op, Args...>::value;
}

#endif // !NA_TYPE_TRAITS_IS_DETECTED_H
