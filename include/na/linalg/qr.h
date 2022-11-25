#ifndef NA_LINALG_QR_H
#define NA_LINALG_QR_H
#include "Eigen/Core"

namespace na
{
	namespace linalg
	{
		// Compute the LQ decomposition of the (n,n) real matrix A
		// such that A = L * Q^'.
		//
		// Input parameters:
		//   A: the (n,n) real matrix to be decomposed
		//
		// Output parameters:
		//   Q: an (n,n) orthogonal matrix
		//   L: an (n,n) lower triangular matrix
		void qr_factorize(
			const Eigen::MatrixXd& A,
			Eigen::MatrixXd& Q,
			Eigen::MatrixXd& L);

		// Compute the solution to the system of linear algebraic equations
		// A * x = y, where the decomposition A = L * Q^' is decomposed by a
		// previous call to qr_factorize().
		//
		// Input parameters:
		//   Q: an (n,n) orthogonal matrix
		//   L: an (n,n) lower triangular matrix
		//   y: the right-hand side of the system of linear algebraic equations
		//
		// Output parameters:
		//   x: the solution to the system of linear algebraic equations
		void qr_solve(
			const Eigen::MatrixXd& Q,
			const Eigen::MatrixXd& L,
			const Eigen::VectorXd& y,
			Eigen::VectorXd& x);

		// Compute the LQ decomposition of the (n,n) complex matrix A
		// such that A = L * Q^*.
		//
		// Input parameters:
		//   A: the (n,n) complex matrix to be decomposed
		//
		// Output parameters:
		//   Q: an (n,n) unitary matrix
		//   L: an (n,n) lower triangular matrix
		void cqr_factorize(
			const Eigen::MatrixXcd& A,
			Eigen::MatrixXcd& Q,
			Eigen::MatrixXcd& L);

		// Compute the solution to the system of linear algebraic equations
		// A * x = y, where the decomposition A = L * Q^* is decomposed by a
		// previous call to cqr_factorize().
		//
		// Input parameters:
		//   Q: an (n,n) unitary matrix
		//   L: an (n,n) lower triangular matrix
		//   y: the right-hand side of the system of linear algebraic equations
		//
		// Output parameters:
		//   x: the solution to the system of linear algebraic equations
		void cqr_solve(
			const Eigen::MatrixXcd& Q,
			const Eigen::MatrixXcd& L,
			const Eigen::VectorXcd& y,
			Eigen::VectorXcd& x);
	}
}

#endif // !NA_LINALG_QR_H
