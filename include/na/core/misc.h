#ifndef NA_CORE_MISC_H
#define NA_CORE_MISC_H
#include <algorithm>

namespace na
{
	template <typename Scalar>
	Scalar clamp(const Scalar& val, const Scalar& min, const Scalar& max)
	{
		return std::max<Scalar>(std::min<Scalar>(val, max), min);
	}

	template <typename Scalar>
	Scalar lerp(const Scalar& a, const Scalar& b, const Scalar& t)
	{
		return a + (b - a) * t;
	}

	template <typename Scalar>
	Scalar bilerp(const Scalar& a, const Scalar& b, const Scalar& c, const Scalar& d, const Scalar& u, const Scalar& v)
	{
		return lerp(lerp(a, b, u), lerp(c, d, u), v);
	}
}

#endif // !NA_CORE_MISC_H
