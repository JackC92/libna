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
		inline Eigen::Quaterniond RotationXY()
		{
			return Eigen::Quaterniond(M_SQRT1_2, 0.0, 0.0, M_SQRT1_2);
		}

		inline Eigen::Quaterniond RotationXZ()
		{
			return Eigen::Quaterniond(M_SQRT1_2, 0.0, -M_SQRT1_2, 0.0);
		}

		inline Eigen::Quaterniond RotationYX()
		{
			return Eigen::Quaterniond(M_SQRT1_2, 0.0, 0.0, -M_SQRT1_2);
		}

		inline Eigen::Quaterniond RotationYZ()
		{
			return Eigen::Quaterniond(M_SQRT1_2, M_SQRT1_2, 0.0, 0.0);
		}

		inline Eigen::Quaterniond RotationZX()
		{
			return Eigen::Quaterniond(M_SQRT1_2, 0.0, M_SQRT1_2, 0.0);
		}

		inline Eigen::Quaterniond RotationZY()
		{
			return Eigen::Quaterniond(M_SQRT1_2, -M_SQRT1_2, 0.0, 0.0);
		}

		inline Eigen::Quaterniond ReflectionX()
		{
			return Eigen::Quaterniond(0.0, 1.0, 0.0, 0.0);
		}

		inline Eigen::Quaterniond ReflectionY()
		{
			return Eigen::Quaterniond(0.0, 0.0, 1.0, 0.0);
		}

		inline Eigen::Quaterniond ReflectionY()
		{
			return Eigen::Quaterniond(0.0, 0.0, 0.0, 1.0);
		}

		// This function is equivalent to Eigen::Quaterniond(Eigen::AngleAxisd(angle_axis.norm(), angle_axis.normalized())).
		// This function is included to have comparable results with Project Chrono.
		template <typename Derived, typename Real = typename Derived::Scalar, int Options = 0>
		inline Eigen::Quaternion<Real, Options> from_rotation_vector(const Eigen::MatrixBase<Derived>& angle_axis)
		{
			static_assert(std::is_floating_point_v<Real>, "from_rotation_vector: Real must be a floating-point type");
			static_assert(Derived::IsVectorAtCompileTime && (Derived::SizeAtCompileTime == 3), "from_rotation_vector: Derived must have compatible size");
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

		template <typename Real, int Options = 0>
		inline Eigen::Quaternion<Real, Options> from_angle_axis(const Real angle, const Eigen::Vector3<Real>& axis)
		{
			static_assert(std::is_floating_point_v<Real>, "from_angle_axis: Real must be a floating-point type");
			Real half_angle = static_cast<Real>(0.5) * angle;
			Real sin_half = std::sin(half_angle);
			// The vector axis is expected to be a unit vector.
			return Eigen::Quaternion<Real, Options>(std::cos(half_angle), axis(0) * sin_half, axis(1) * sin_half, axis(2) * sin_half);
		}

		template <typename Real, int Options = 0>
		inline Eigen::Quaternion<Real, Options> from_angleX(const Real angleX)
		{
			static_assert(std::is_floating_point_v<Real>, "from_angleX: Real must be a floating-point type");
			Real half_angle = static_cast<Real>(0.5) * angleX;
			return Eigen::Quaternion<Real, Options>(std::cos(half_angle), std::sin(half_angle), 0.0, 0.0);
		}

		template <typename Real, int Options = 0>
		inline Eigen::Quaternion<Real, Options> from_angleY(const Real angleY)
		{
			static_assert(std::is_floating_point_v<Real>, "from_angleY: Real must be a floating-point type");
			Real half_angle = static_cast<Real>(0.5) * angleY;
			return Eigen::Quaternion<Real, Options>(std::cos(half_angle), std::sin(half_angle), 0.0, 0.0);
		}

		template <typename Real, int Options = 0>
		inline Eigen::Quaternion<Real, Options> from_angleZ(const Real angleZ)
		{
			static_assert(std::is_floating_point_v<Real>, "from_angleZ: Real must be a floating-point type");
			Real half_angle = static_cast<Real>(0.5) * angleZ;
			return Eigen::Quaternion<Real, Options>(std::cos(half_angle), std::sin(half_angle), 0.0, 0.0);
		}

		/*
		* The following functions are useful for multi-body physics when considering quaternion-valued functions defined over a time interval.
		* Notice that addition, subtraction and scaling of quaternions are not available in Eigen because Eigen::Quaternion was only intended to represent 3D rotations, see https://eigen.tuxfamily.org/bz/show_bug.cgi?id=560.
		*/

		// Compute the quaternion first derivative from the vector of angular speed given in absolute coordinates.
		template <typename Real, int Options>
		inline Eigen::Quaternion<Real, Options> Qdt_from_Wabs(const Eigen::Vector3<Real>& w, const Eigen::Quaternion<Real, Options>& q)
		{
			return Eigen::Quaternion<Real, Options>(0.0, w(0) * static_cast<Real>(0.5), w(1) * static_cast<Real>(0.5), w(2) * static_cast<Real>(0.5)) * q;
		}

		// Compute the quaternion first derivative from the vector of angular speed given in local coordinates.
		template <typename Real, int Options>
		inline Eigen::Quaternion<Real, Options> Qdt_from_Wloc(const Eigen::Vector3<Real>& w, const Eigen::Quaternion<Real, Options>& q)
		{
			return q * Eigen::Quaternion<Real, Options>(0.0, w(0) * static_cast<Real>(0.5), w(1) * static_cast<Real>(0.5), w(2) * static_cast<Real>(0.5));
		}

		// Compute the quaternion second derivative from the vector of angular speed given in absolute coordinates.
		template <typename Real, int Options>
		inline Eigen::Quaternion<Real, Options> Qdtdt_from_Aabs(
			const Eigen::Vector3<Real>& a,
			const Eigen::Quaternion<Real, Options>& q,
			const Eigen::Quaternion<Real, Options>& q_dt)
		{
			return Eigen::Quaternion<Real, Options>(0.0, a(0) * static_cast<Real>(0.5), a(1) * static_cast<Real>(0.5), a(2) * static_cast<Real>(0.5)) * q + q_dt * q.conjugate() * q_dt;
		}

		// Compute the quaternion second derivative from the vector of angular speed given in local coordinates.
		template <typename Real, int Options>
		inline Eigen::Quaternion<Real, Options> Qdtdt_from_Aloc(
			const Eigen::Vector3<Real>& a,
			const Eigen::Quaternion<Real, Options>& q,
			const Eigen::Quaternion<Real, Options>& q_dt)
		{
			return q * Eigen::Quaternion<Real, Options>(0.0, a(0) * static_cast<Real>(0.5), a(1) * static_cast<Real>(0.5), a(2) * static_cast<Real>(0.5)) + q_dt * q.conjugate() * q_dt;
		}
	}
}

#endif // !NA_GEOMETRY_QUATERNION_H
