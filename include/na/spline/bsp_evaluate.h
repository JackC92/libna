#ifndef NA_SPLINE_BSP_EVALUATE_H
#define NA_SPLINE_BSP_EVALUATE_H
#include "Eigen/Core"
#include "Eigen/Sparse"

namespace na
{
	namespace bspline
	{
		Eigen::VectorXd evaluate(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index degree,
			const Eigen::Index m,
			const double u);
		
		Eigen::VectorXd evaluate_basis(
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index degree,
			const Eigen::Index m,
			const Eigen::Index mu,
			const double u);

		Eigen::SparseMatrix<double, Eigen::RowMajor> evaluate_basis(
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index degree,
			const Eigen::Index m,
			const Eigen::Ref<const Eigen::ArrayXi>& mu,
			const Eigen::Ref<const Eigen::ArrayXd>& u);

		double length(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index degree,
			const double u0,
			const double u1,
			const double tol = 1e-12);
		
		// Compute the curvature binormal of the B-spline curve in 3-space.
		Eigen::Vector3d curvature_binormal(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index degree,
			const double u);
		
		// Compute the curvature binormal scaled by the norm of the tangent
		// This is usually used for solving the parallel transport ODE parametrized by the arc-length parameter.
		Eigen::Vector3d scaled_curvature_binormal(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index degree,
			const double u);
	}
}

#endif // !NA_SPLINE_BSP_EVALUATE_H
