#ifndef NA_LINALG_HESSENBERG_H
#define NA_LINALG_HESSENBERG_H
#include <complex>
#include "Eigen/Core"

namespace na
{
	namespace linalg
	{
		// Compute the eigenvalues of the lower Hessenberg matrix B of
			// the form
			//              B = A + P \circ Q^*,                                   (1)
			// where A is a Hermitian matrix, and P, Q are a pair of vectors.
			// Note that this is a very specialized subroutine, in that the
			// matrix B is BOTH of the form (1) and lower Hessenber, so that
			//              B(i, j) = 0                                            (2)
			// for all
			//              j > i + 1.                                             (3)
			// Note also that, for matrices of this form, A is determined
			// entirely by its diagonal and superdiagonal together with the
			// vectors P and Q.
		void herm_p_rank1(
			Eigen::VectorXcd& diag,
			Eigen::VectorXcd& supdiag,
			Eigen::VectorXcd& p,
			Eigen::VectorXcd& q,
			const Eigen::Index n,
			const double eps);
	}
}

#endif // !NA_LINALG_HESSENBERG_H
