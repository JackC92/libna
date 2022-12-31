#ifndef NA_CORE_STRING_H
#define NA_CORE_STRING_H
#include <cstdio>
#include <string>

namespace na
{
	template <typename... Args>
	std::string format(const std::string& fmt, Args... args)
	{
		int size = std::snprintf(nullptr, 0, fmt.c_str(), args...) + 1;
		std::string output(size + 1, '\0');
		std::sprintf(&output[0], fmt.c_str(), args...);
		return output;
	}
}

#endif // !NA_CORE_STRING_H
