#ifndef NA_GEOMETRY_FRAME_MAPPING_H
#define NA_GEOMETRY_FRAME_MAPPING_H
#include "Eigen/Core"

namespace na
{
	// This function is based on the paper
	//   C. Meier, A. Popp, and W. A. Wall, "An objective 3D large deformation finite element formulation for geometrically exact curved Kirchhoff rods",
	//   Computer methods in applied mechanics and engineering, 2014
	// This is the same as discrete parallel transport, and is equivalent to Eigen::Quaterniond::FromTwoVectors(t0, t1)._transformVector(v).
	Eigen::Vector3d smallest_rotation(
		const Eigen::Vector3d& t0,
		const Eigen::Vector3d& t1,
		const Eigen::Vector3d& v);

	// This function is based on the paper
	//   W. Wang, B. Juttler, D. Zheng, and Y. Liu, "Computation of Rotation Minimizing Frames",
	//   ACM Transactions on graphics, Vol. 27, No. 1, Article 2, 2008
	Eigen::Vector3d double_reflection(
		const Eigen::Vector3d& x0,
		const Eigen::Vector3d& t0,
		const Eigen::Vector3d& x1,
		const Eigen::Vector3d& t1,
		const Eigen::Vector3d& v);
}

#endif // !NA_GEOMETRY_FRAME_MAPPING_H
