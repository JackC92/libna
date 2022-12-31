#include "na/core/matrix.h"
#include "Eigen/Core"

namespace na
{
	Eigen::Matrix3d hat_matrix(const Eigen::Vector3d& vec)
	{
		Eigen::Matrix3d mat;
		mat << 0.0, -vec(2), vec(1),
			vec(2), 0.0, -vec(0),
			-vec(1), vec(0), 0.0;
		return mat;
	}
}
