#ifndef NA_CHEBYSHEV_H
#define NA_CHEBYSHEV_H
#include "Eigen/Core"

namespace na
{
	namespace chebyshev
	{
		void quadrature(
			const Eigen::Index n,
			Eigen::VectorXd& xscheb,
			Eigen::VectorXd& whtscheb);

		Eigen::VectorXd expansion_coefficients(
			const Eigen::Index n,
			const Eigen::VectorXd& vals);

		double evaluate_expansion(
			const Eigen::Index n,
			const Eigen::VectorXd& coefs,
			const double x);

		void coefficient_matrix(
			const Eigen::Index n,
			const Eigen::VectorXd& xscheb,
			const Eigen::VectorXd& whtscheb,
			Eigen::MatrixXd& umatr,
			Eigen::MatrixXd& vmatr);
	}

	namespace chebyshev_pract
	{
		void quadrature(
			const Eigen::Index n,
			Eigen::VectorXd& xscheb,
			Eigen::VectorXd& whtscheb);

		void coefficient_matrix(
			const Eigen::Index n,
			const Eigen::VectorXd& xscheb,
			const Eigen::VectorXd& whtscheb,
			Eigen::MatrixXd& umatr,
			Eigen::MatrixXd& vmatr);

		Eigen::MatrixXd differentiation_matrix(
			const Eigen::Index n,
			const Eigen::VectorXd& xscheb);
	}
}

#endif // !NA_CHEBYSHEV_H
