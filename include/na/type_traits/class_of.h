#ifndef NA_TYPE_TRAITS_CLASS_OF_H
#define NA_TYPE_TRAITS_CLASS_OF_H

namespace na
{
	template <typename T>
	struct class_of {};

	template <typename ReturnType, typename Class>
	struct class_of<ReturnType(Class::*)>
	{
		using type = Class;
	};
	
	template <typename T>
	using class_of_t = typename class_of<T>::type;
}

#endif // !NA_TYPE_TRAITS_CLASS_OF_H
