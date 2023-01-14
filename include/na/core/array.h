#ifndef NA_CORE_ARRAY_H
#define NA_CORE_ARRAY_H
#include <algorithm>
#include "Eigen/Core"

namespace na
{
	template <typename Scalar>
	inline void deduplicate(
		Eigen::ArrayX<Scalar>& arr,
		const bool sort = false)
	{
		if (arr.size() == 0)
		{
			return;
		}
		if (sort)
		{
			std::sort<Scalar*>(arr.data(), arr.data() + arr.size());
		}
		arr.conservativeResize(std::unique<Scalar*>(arr.data(), arr.data() + arr.size()) - arr.data());
	}

	template <typename Derived>
	inline void repelem(
		const Eigen::ArrayBase<Derived>& arr,
		const Eigen::Index n,
		Eigen::ArrayX<typename Derived::Scalar>& res)
	{
		static_assert(Derived::IsVectorAtCompileTime, "repelem: Derived must be a vector at compile time");
		res = arr.reshaped().transpose().replicate(n, 1).reshaped();
	}

	template <typename Derived>
	inline Eigen::ArrayX<typename Derived::Scalar> repelem(
		const Eigen::ArrayBase<Derived>& vec,
		const Eigen::Index n)
	{
		Eigen::ArrayX<typename Derived::Scalar> res;
		repelem(vec, n, res);
		return res;
	}
}

#endif // !NA_CORE_ARRAY_H
