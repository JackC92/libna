#include "na/bspline.h"
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
			const double s)
		{
			Eigen::VectorXd b;
			int mu;
			eval_basis(T, deg, m, s, mu, b);
			return q.middleRows(mu - deg, deg + 1).transpose() * b;
		}

		void eval_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s,
			int& mu,
			Eigen::VectorXd& b)
		{
			// n is the number of control points needed to be compatible with T.
			const int n = T.size() - deg - 1;

			assert((n >= deg + 1) && "T must contain at least 2 * deg + 2 knots");
			assert(((T(deg) <= s) && (s <= T(n))) && "s must be between [T(deg), T(n)]");
			assert((0 <= m) && "m must be a non-negative integer");

			// The index of the left endpoint of the interval [T(mu), T(mu) + 1) containing s.
			// The only exception is s = T(n): mu = n - 1 in this case, and the basis functions can be computed correctly using mu = n - 1.
			mu = deg;
			while (mu < n - 1)
			{
				if ((T(mu) <= s) && (s < T(mu + 1)))
				{
					break;
				}
				mu += 1;
			}

			// s is in the support of exactly deg + 1 B-spline basis functions.
			b.setZero(deg + 1);
			if (m <= deg)
			{
				b(deg) = 1.0;

				double w1, w2;
				// B_{0} = 1 is a 1x1 matrix and B_{r}(s) = B_{r - 1}(s) * R_{r}(s).
				// Post-multiply by the matrix R_{r}(s).
				for (int r = 1; r <= deg - m; ++r)
				{
					int k = mu - r + 1;
					w2 = (T(k + r) - s) / (T(k + r) - T(k));
					b(deg - r) = w2 * b(deg - r + 1);
					for (int i = deg - r + 1; i <= deg - 1; ++i)
					{
						k += 1;
						w1 = w2;
						w2 = (T(k + r) - s) / (T(k + r) - T(k));
						b(i) = (1.0 - w1) * b(i) + w2 * b(i + 1);
					}
					b(deg) *= 1.0 - w2;
				}

				// Post-multiply by the matrix r * DR_r, which is independent of s because the entries of R_r(s) are all linear in s.
				for (int r = deg - m + 1; r <= deg; ++r)
				{
					int k = mu - r + 1;
					b(deg - r) = -b(deg - r + 1) / (T(k + r) - T(k));
					for (int i = deg - r + 1; i <= deg - 1; ++i)
					{
						b(i) = b(i) / (T(k + r) - T(k)) - b(i + 1) / (T(k + r + 1) - T(k + 1));
						k += 1;
					}
					b(deg) /= T(k + r) - T(k);
					b *= r;
				}
			}
		}

		Eigen::VectorXd eval_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s)
		{
			Eigen::VectorXd b;
			int mu;
			eval_basis(T, deg, m, s, mu, b);
			return b;
		}

		Eigen::MatrixXd eval_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const Eigen::VectorXd& s)
		{
			// n is the number of control points needed to be compatible with T.
			const int n = T.size() - deg - 1;

			// b.row(i) contains the value of all n basis function evald at s(i).
			const int num_pts = s.size();
			Eigen::MatrixXd b(num_pts, n);
			b.setZero();

			int mu;
			Eigen::VectorXd temp;
			for (int i = 0; i < num_pts; ++i)
			{
				eval_basis(T, deg, m, s(i), mu, temp);
				b.block(i, mu - deg, 1, deg + 1) = temp.transpose();
			}
			return b;
		}
	}
}
