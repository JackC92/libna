#ifndef NA_SPLINE_BSP_EVALUATE_H
#define NA_SPLINE_BSP_EVALUATE_H
#include "Eigen/Core"

namespace na
{
	namespace bspline
	{
		/*
		* Note that Eigen::Ref can avoid data copying only if the input argument has the same storage order as the template argument for Eigen::Ref.
		* For example, Eigen::Transpose<Eigen::MatrixXd> will result in data copying into a temporary object because of the different storage order.
		*/

		Eigen::VectorXd evaluate(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::VectorXd>& knots,
			const Eigen::Index degree,
			const Eigen::Index m,
			const double u);

		Eigen::VectorXd evaluate_basis(
			const Eigen::Ref<const Eigen::VectorXd>& knots,
			const Eigen::Index degree,
			const Eigen::Index m,
			const double u);

		Eigen::MatrixXd evaluate_basis(
			const Eigen::Ref<const Eigen::VectorXd>& knots,
			const Eigen::Index degree,
			const Eigen::Index m,
			const Eigen::Ref<const Eigen::VectorXd>& u);

		double length(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::VectorXd>& knots,
			const Eigen::Index degree,
			const double u0,
			const double u1,
			const double tol = 1e-6);

		// Compute the curvature binormal of the B-spline curve described by the input parameters.
		Eigen::Vector3d curvature_binormal(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::VectorXd>& knots,
			const Eigen::Index degree,
			const double u);

		// Return the curvature binormal scaled by the norm of the tangent.
		// This is usually used for solving the parallel transport ODE.
		Eigen::Vector3d scaled_curvature_binormal(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::VectorXd>& knots,
			const Eigen::Index degree,
			const double u);
	}
}

#endif // !NA_SPLINE_BSP_EVALUATE_H
