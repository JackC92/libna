#ifndef NA_SPLINE_BZR_EVALUATE_H
#define NA_SPLINE_BZR_EVALUATE_H
#include "Eigen/Core"

namespace na
{
	namespace bezier
	{
		Eigen::Vector3d evaluate(
			const Eigen::Ref<const Eigen::MatrixX3d>& coefs,
			const double u);

		Eigen::Vector3d evaluate(
			const Eigen::Ref<const Eigen::MatrixX4d>& coefs,
			const double u);
	}
}

#endif // !NA_SPLINE_BZR_EVALUATE_H
