#ifndef NA_BSPLINE_H
#define NA_BSPLINE_H
#include "Eigen/Core"

namespace na
{
	namespace bspline
	{
		Eigen::Vector3d evaluate(
			const Eigen::MatrixXd& q,
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s);

		void evaluate_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s,
			int& mu,
			Eigen::VectorXd& b);

		Eigen::VectorXd evaluate_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s);

		Eigen::MatrixXd evaluate_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const Eigen::VectorXd& s);
	}
}

#endif // !NA_BSPLINE_H
