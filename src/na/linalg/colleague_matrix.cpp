#include "na/linalg/colleague_matrix.h"
#include <cassert>
#include <cmath>
#include <complex>
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "na/special/chebyshev.h"
#include "na/linalg/hessenberg.h"

namespace na
{
	namespace linalg
	{
		Eigen::SparseMatrix<double> colleague_matrix(
			const Eigen::Index n,
			const Eigen::VectorXd& coefs)
		{
			assert((n >= 1) && "colleague_matrix: n must be a positive integer");
			assert((coefs.size() == n + 1) && "colleague_matrix: n must be the order of the Chebyshev expansion");

			Eigen::SparseMatrix<double> cmatr(n, n);
			cmatr.reserve(Eigen::VectorXi::Constant(n, 3));
			cmatr.insert(1, 0) = std::sqrt(0.5);
			cmatr.insert(0, 1) = std::sqrt(0.5);
			for (Eigen::Index i = 1; i < n - 1; ++i)
			{
				cmatr.insert(i, i + 1) = 0.5;
				cmatr.insert(i + 1, i) = 0.5;
			}
			for (Eigen::Index i = 0; i < n; ++i)
			{
				cmatr.coeffRef(n - 1, i) += -coefs(i) / coefs(n) * (i == 0 ? std::sqrt(0.5) : 0.5);
			}
			return cmatr;
		}

		Eigen::SparseMatrix<Eigen::dcomplex> colleague_matrix(
			const Eigen::Index n,
			const Eigen::VectorXcd& coefs)
		{
			assert((n >= 1) && "colleague_matrix: n must be a positive integer");
			assert((coefs.size() == n + 1) && "colleague_matrix: n must be the order of the Chebyshev expansion");

			Eigen::SparseMatrix<Eigen::dcomplex> cmatr(n, n);
			cmatr.reserve(Eigen::VectorXi::Constant(n, 3));
			cmatr.insert(1, 0) = std::sqrt(0.5);
			cmatr.insert(0, 1) = std::sqrt(0.5);
			for (Eigen::Index i = 1; i < n - 1; ++i)
			{
				cmatr.insert(i, i + 1) = 0.5;
				cmatr.insert(i + 1, i) = 0.5;
			}
			for (Eigen::Index i = 0; i < n; ++i)
			{
				cmatr.coeffRef(n - 1, i) += -std::conj(coefs(i) / coefs(n)) * (i == 0 ? std::sqrt(0.5) : 0.5);
			}
			return cmatr;
		}

		Eigen::VectorXcd colleague_roots(
			const Eigen::Index n,
			const Eigen::VectorXcd& coefs,
			const double eps)
		{
			assert((n >= 1) && "colleague_roots: n must be a positive integer");

			if (n == 1)
			{
				return Eigen::VectorXcd::Constant(1, -coefs(0) / coefs(1));
			}

			Eigen::VectorXcd roots = Eigen::VectorXcd::Zero(n);
			Eigen::VectorXcd supdiag = Eigen::VectorXcd::Constant(n - 1, 0.5);
			Eigen::VectorXcd p = Eigen::VectorXcd::Unit(n, n - 1);
			Eigen::VectorXcd q = (-0.5 * coefs.head(n) / coefs(n)).conjugate();
			supdiag(0) = std::sqrt(0.5);
			q(0) *= std::sqrt(2.0);
			na::linalg::herm_p_rank1(roots, supdiag, p, q, n, eps);
			return roots;
		}

		Eigen::VectorXd colleague_roots_m1p1(
			const Eigen::Index n,
			const Eigen::VectorXd& coefs,
			const double eps,
			const double delta,
			const double coff)
		{
			assert((n >= 1) && "colleague_roots_m1p1: n must be a positive integer");

			Eigen::VectorXd roots(n);
			Eigen::VectorXcd croots = colleague_roots(n, coefs, eps);
			// Select all the complex roots in the rectangle [-1-delta, 1+delta] x [-delta, delta] \subset \C
			int count = 0;
			for (Eigen::Index i = 0; i < n; ++i)
			{
				if ((std::abs(croots(i).real()) > 1.0 + delta) || (std::abs(croots(i).imag()) > delta))
				{
					continue;
				}
				roots(count) = croots(i).real();
				count += 1;
			}
			// Discard any roots with a large error
			double coefsize = coefs.norm();
			int nroots = 0;
			for (Eigen::Index i = 0; i < count; ++i)
			{
				if (std::abs(na::chebyshev::evaluate_expansion(n, coefs, roots(i))) > coff * coefsize)
				{
					continue;
				}
				roots(nroots) = roots(i);
				nroots += 1;
			}
			return roots.head(nroots);
		}
	}
}
