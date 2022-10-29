#include "na/geometry/frame_mapping.h"
#include "Eigen/Core"

namespace na
{
	void make_basis(
		const Eigen::Vector3d& D3,
		Eigen::Ref<Eigen::Vector3d> D1,
		Eigen::Ref<Eigen::Vector3d> D2)
	{
		if (D3(1) * D3(1) + D3(2) * D3(2) > 1e-12)
		{
			D1 = Eigen::Vector3d(0.0, D3(2), -D3(1)).normalized();
		}
		else
		{
			D1 = Eigen::Vector3d(-D3(2), 0.0, D3(0)).normalized();
		}
		D2 = D3.cross(D1);
	}

	Eigen::Vector3d smallest_rotation(
		const Eigen::Vector3d& t0,
		const Eigen::Vector3d& t1,
		const Eigen::Vector3d& v)
	{
		// t0 and t1 must be unit tangent vectors.
		return v - v.dot(t1) / (1.0 + t0.dot(t1)) * (t0 + t1);
	}

	Eigen::Vector3d double_reflection(
		const Eigen::Vector3d& x0,
		const Eigen::Vector3d& t0,
		const Eigen::Vector3d& x1,
		const Eigen::Vector3d& t1,
		const Eigen::Vector3d& v)
	{
		// t0 and t1 must be unit tangent vectors.
		Eigen::Vector3d v1 = x1 - x0;
		Eigen::Vector3d rL = v - 2.0 / v1.squaredNorm() * v1.dot(v) * v1;
		Eigen::Vector3d tL = t0 - 2.0 / v1.squaredNorm() * v1.dot(t0) * v1;
		Eigen::Vector3d v2 = t1 - tL;
		return rL - 2.0 / v2.squaredNorm() * v2.dot(rL) * v2;
	}
}
