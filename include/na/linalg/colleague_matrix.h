#ifndef NA_LINALG_COLLEAGUE_MATRIX_H
#define NA_LINALG_COLLEAGUE_MATRIX_H
#include "Eigen/Core"
#include "Eigen/Sparse"

namespace na
{
	namespace linalg
	{
		Eigen::SparseMatrix<double> colleague_matrix(
			const Eigen::Index n,
			const Eigen::VectorXd& coefs);

		Eigen::SparseMatrix<Eigen::dcomplex> colleague_matrix(
			const Eigen::Index n,
			const Eigen::VectorXcd& coefs);

		Eigen::VectorXcd colleague_roots(
			const Eigen::Index n,
			const Eigen::VectorXcd& coefs,
			const double eps);

		Eigen::VectorXd colleague_roots_m1p1(
			const Eigen::Index n,
			const Eigen::VectorXd& coefs,
			const double eps,
			const double delta,
			const double coff);
	}
}

#endif // !NA_LINALG_COLLEAGUE_MATRIX_H
