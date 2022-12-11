#ifndef NA_CORE_VECTOR_H
#define NA_CORE_VECTOR_H
#include <algorithm>
#include "Eigen/Core"

namespace na
{
	template <typename Scalar>
	inline void deduplicate(
		Eigen::VectorX<Scalar>& vec,
		const bool sort = false)
	{
		if (vec.size() == 0)
		{
			return;
		}
		if (sort)
		{
			std::sort<Scalar*>(vec.data(), vec.data() + vec.size());
		}
		vec.conservativeResize(std::unique<Scalar*>(vec.data(), vec.data() + vec.size()) - vec.data());
	}

	template <typename Derived>
	inline void repelem(
		const Eigen::MatrixBase<Derived>& vec,
		const Eigen::Index n,
		Eigen::VectorX<typename Derived::Scalar>& res)
	{
		static_assert(Derived::IsVectorAtCompileTime, "repelem: Derived must be a vector at compile time");
		res = vec.reshaped().transpose().replicate(n, 1).reshaped();
	}
	
	template <typename Derived>
	inline Eigen::VectorX<typename Derived::Scalar> repelem(
		const Eigen::MatrixBase<Derived>& vec,
		const Eigen::Index n)
	{
		Eigen::VectorX<typename Derived::Scalar> res;
		repelem(vec, n, res);
		return res;
	}
}

#endif // !NA_CORE_VECTOR_H
