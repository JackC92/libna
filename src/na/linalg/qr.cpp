#include "na/linalg/qr.h"
#include <cassert>
#include <cmath>
#include <complex>
#include "Eigen/Core"

namespace na
{
	namespace internal
	{
		namespace linalg
		{
			// Compute a (2,2) orthogonal matrix u which satisfies
			//         u * (x, y)^' = (z, 0)^'
			inline Eigen::Matrix2d qr_rotation_matrix(
				const double x,
				const double y,
				const double size)
			{
				Eigen::Matrix2d u;
				double d = x * x + y * y;
				if (d < size * 1e-66)
				{
					u = Eigen::Vector2d::Ones().asDiagonal();
				}
				else
				{
					d = std::sqrt(d);
					u(0, 0) = -x / d;
					u(1, 0) = -y / d;
					u(0, 1) = u(1, 0);
					u(1, 1) = -u(0, 0);
				}
				return u;
			}

			// Compute a (2,2) unitary matrix u which satisfies
			//         u * (x, y)^' = (z, 0)^'
			inline Eigen::Matrix2cd cqr_rotation_matrix(
				const Eigen::dcomplex& x,
				const Eigen::dcomplex& y,
				const double size)
			{
				Eigen::Matrix2cd u;
				double d = std::norm(x) + std::norm(y);
				if (d < size * 1e-66)
				{
					u = Eigen::Vector2cd::Ones().asDiagonal();
				}
				else
				{
					d = std::sqrt(d);
					u(0, 1) = -y / d;
					u(1, 1) = x / d;
					u(0, 0) = -std::conj(u(1, 1));
					u(1, 0) = std::conj(u(0, 1));
				}
				return u;
			}
			
			// Apply the inverse of L to y, and multiply the result by Q on the left.
			template <typename Scalar>
			inline void qr_apply(
				const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& Q,
				const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& L,
				const Eigen::Vector<Scalar, Eigen::Dynamic>& y,
				Eigen::Vector<Scalar, Eigen::Dynamic>& x)
			{
				const int n = y.size();
				x = y;
				x(0) /= L(0, 0);
				for (int i = 1; i < n; ++i)
				{
					x.segment(i, n - i) -= x(i - 1) * L.col(i - 1).segment(i, n - i);
					x(i) /= L(i, i);
				}
				x = Q * x;
			}
		}
	}

	namespace linalg
	{
		void qr_factorize(
			const Eigen::MatrixXd& A,
			Eigen::MatrixXd& Q,
			Eigen::MatrixXd& L)
		{
			assert((A.rows() == A.cols()) && "A must be a square matrix");

			const int n = A.rows();
			const double size = A.squaredNorm();
			L = A;
			Q = Eigen::VectorXd::Ones(n).asDiagonal();

			for (int i = 0; i < n - 1; ++i)
			{
				for (int j = n - 1; j > i; --j)
				{
					Eigen::Matrix2d u = na::internal::linalg::qr_rotation_matrix(L(i, i), L(i, j), size);
					L(Eigen::seq(i, n - 1), { i, j }) = L(Eigen::seq(i, n - 1), { i, j }) * u;
					Q(Eigen::placeholders::all, { i, j }) = Q(Eigen::placeholders::all, { i, j }) * u;
				}
			}
		}

		void qr_solve(
			const Eigen::MatrixXd& Q,
			const Eigen::MatrixXd& L,
			const Eigen::VectorXd& y,
			Eigen::VectorXd& x)
		{
			na::internal::linalg::qr_apply(Q, L, y, x);
		}

		void cqr_factorize(
			const Eigen::MatrixXcd& A,
			Eigen::MatrixXcd& Q,
			Eigen::MatrixXcd& L)
		{
			assert((A.rows() == A.cols()) && "A must be a square matrix");

			const int n = A.rows();
			const double size = A.squaredNorm();
			L = A;
			Q = Eigen::VectorXcd::Ones(n).asDiagonal();

			for (int i = 0; i < n - 1; ++i)
			{
				for (int j = n - 1; j > i; --j)
				{
					Eigen::Matrix2cd u = na::internal::linalg::cqr_rotation_matrix(L(i, i), L(i, j), size);
					L(Eigen::seq(i, n - 1), { i, j }) = L(Eigen::seq(i, n - 1), { i, j }) * u;
					Q(Eigen::placeholders::all, { i, j }) = Q(Eigen::placeholders::all, { i, j }) * u;
				}
			}
		}

		void cqr_solve(
			const Eigen::MatrixXcd& Q,
			const Eigen::MatrixXcd& L,
			const Eigen::VectorXcd& y,
			Eigen::VectorXcd& x)
		{
			na::internal::linalg::qr_apply(Q, L, y, x);
		}
	}
}