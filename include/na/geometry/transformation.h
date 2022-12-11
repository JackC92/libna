#ifndef NA_GEOMETRY_TRANSFORMATION_H
#define NA_GEOMETRY_TRANSFORMATION_H
#include <type_traits>
#include "Eigen/Core"
#include "Eigen/Geometry"

namespace na
{
	template <typename Scalar, int Dim>
	Eigen::Matrix<Scalar, Dim + 1, Dim + 1> homogeneous_translation(const Eigen::Vector<Scalar, Dim>& vec)
	{
		static_assert(std::is_floating_point_v<Scalar>, "homogeneous_translation: Scalar must be a floating-point type");
		Eigen::Matrix<Scalar, Dim + 1, Dim + 1> mat = Eigen::Matrix<Scalar, Dim + 1, Dim + 1>::Identity();
		mat.template block<Dim, 1>(0, Dim) = vec;
		return mat;
	}

	template <typename Scalar, int Dim>
	Eigen::Matrix<Scalar, Dim + 1, Dim + 1> homogeneous_scaling(const Eigen::Vector<Scalar, Dim>& vec)
	{
		static_assert(std::is_floating_point_v<Scalar>, "homogeneous_scaling: Scalar must be a floating-point type");
		Eigen::Matrix<Scalar, Dim + 1, Dim + 1> mat = Eigen::Matrix<Scalar, Dim + 1, Dim + 1>::Identity();
		mat.template block<Dim, Dim>(0, 0).diagonal() = vec;
		return mat;
	}

	template <typename Scalar>
	Eigen::Matrix3<Scalar> homogeneous_rotation(const Scalar& angle)
	{
		static_assert(std::is_floating_point_v<Scalar>, "homogeneous_rotation: Scalar must be a floating-point type");
		Eigen::Matrix3<Scalar> mat = Eigen::Matrix3<Scalar>::Identity();
		mat.template block<2, 2>(0, 0) = Eigen::Rotation2D<Scalar>(angle).toRotationMatrix();
		return mat;
	}

	template <typename Scalar>
	Eigen::Matrix4<Scalar> homogeneous_rotation(const Scalar& angle, const Eigen::Vector3<Scalar>& axis)
	{
		static_assert(std::is_floating_point_v<Scalar>, "homogeneous_rotation: Scalar must be a floating-point type");
		Eigen::Matrix4<Scalar> mat = Eigen::Matrix4<Scalar>::Identity();
		mat.template block<3, 3>(0, 0) = Eigen::AngleAxis<Scalar>(angle, axis).toRotationMatrix();
		return mat;
	}

	template <typename Scalar>
	Eigen::Matrix4<Scalar> homogeneous_rotation(const Eigen::Quaternion<Scalar>& quat)
	{
		static_assert(std::is_floating_point_v<Scalar>, "homogeneous_rotation: Scalar must be a floating-point type");
		Eigen::Matrix4<Scalar> mat = Eigen::Matrix4<Scalar>::Identity();
		mat.template block<3, 3>(0, 0) = quat.toRotationMatrix();
		return mat;
	}

	template <typename Scalar>
	Eigen::Matrix4<Scalar> homogeneous_rotation_angleX(const Scalar& angle)
	{
		static_assert(std::is_floating_point_v<Scalar>, "homogeneous_rotation_angleX: Scalar must be a floating-point type");
		return homogeneous_rotation(angle, Eigen::Vector3<Scalar>(1.0, 0.0, 0.0));
	}

	template <typename Scalar>
	Eigen::Matrix4<Scalar> homogeneous_rotation_angleY(const Scalar& angle)
	{
		static_assert(std::is_floating_point_v<Scalar>, "homogeneous_rotation_angleY: Scalar must be a floating-point type");
		return homogeneous_rotation(angle, Eigen::Vector3<Scalar>(0.0, 1.0, 0.0));
	}

	template <typename Scalar>
	Eigen::Matrix4<Scalar> homogeneous_rotation_angleZ(const Scalar& angle)
	{
		static_assert(std::is_floating_point_v<Scalar>, "homogeneous_rotation_angleZ: Scalar must be a floating-point type");
		return homogeneous_rotation(angle, Eigen::Vector3<Scalar>(0.0, 0.0, 1.0));
	}
}

#endif // !NA_GEOMETRY_TRANSFORMATION_H
