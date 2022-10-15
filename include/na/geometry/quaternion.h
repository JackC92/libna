#ifndef NA_GEOMETRY_QUATERNION_H
#define NA_GEOMETRY_QUATERNION_H
#include <cmath>
#include <type_traits>
#include "Eigen/Core"
#include "Eigen/Geometry"

template <typename Scalar, int Options>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR const Eigen::Quaternion<Scalar, Options> operator+(
	const Eigen::Quaternion<Scalar, Options>& lhs,
	const Eigen::Quaternion<Scalar, Options>& rhs)
{
	return Eigen::Quaternion<Scalar, Options>(lhs.coeffs() + rhs.coeffs());
}

template <typename Scalar, int Options>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE EIGEN_CONSTEXPR const Eigen::Quaternion<Scalar, Options> operator-(
	const Eigen::Quaternion<Scalar, Options>& lhs,
	const Eigen::Quaternion<Scalar, Options>& rhs)
{
	return Eigen::Quaternion<Scalar, Options>(lhs.coeffs() - rhs.coeffs());
}

namespace na
{
	namespace quaternion
	{
		// This function is equivalent to Eigen::Quaterniond(Eigen::AngleAxisd(angle_axis.norm(), angle_axis.normalized())).
		// This function is included to have comparable results with Project Chrono.
		template <typename Derived, typename Real = typename Derived::Scalar, int Options = 0>
		inline Eigen::Quaternion<Real, Options> from_rotation_vector(const Eigen::MatrixBase<Derived>& angle_axis)
		{
			static_assert(std::is_floating_point_v<Real>, "from_rotation_vector only accepts floating-point types");
			static_assert(Derived::IsVectorAtCompileTime && (Derived::SizeAtCompileTime == 3), "from_rotation_vector only accepts vectors of a specific size");
			Real theta_squared = angle_axis.squaredNorm();
			if (theta_squared > 1e-30)
			{
				Real theta = std::sqrt(theta_squared);
				Real k = std::sin(static_cast<Real>(0.5) * theta) / theta;
				return Eigen::Quaternion<Real, Options>(std::cos(static_cast<Real>(0.5) * theta), angle_axis(0) * k, angle_axis(1) * k, angle_axis(2) * k);
			}
			else
			{
				return Eigen::Quaternion<Real, Options>(static_cast<Real>(1.0), angle_axis(0) * static_cast<Real>(0.5), angle_axis(1) * static_cast<Real>(0.5), angle_axis(2) * static_cast<Real>(0.5));
			}
		}
	}
}

#endif // !NA_GEOMETRY_QUATERNION_H
