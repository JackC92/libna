#ifndef NA_BSPLINE_H
#define NA_BSPLINE_H
#include "Eigen/Core"

namespace na
{
	namespace bspline
	{
		Eigen::Vector3d eval(
			const Eigen::MatrixXd& q,
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s);

		void eval_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s,
			int& mu,
			Eigen::VectorXd& b);

		Eigen::VectorXd eval_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s);

		Eigen::MatrixXd eval_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const Eigen::VectorXd& s);
	}
}

#endif // !NA_BSPLINE_H
