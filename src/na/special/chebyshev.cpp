#include "na/special/chebyshev.h"
#include <cassert>
#include <cmath>
#include "Eigen/Core"

namespace na
{
	namespace chebyshev
	{
		void quadrature(
			const Eigen::Index n,
			Eigen::VectorXd& xscheb,
			Eigen::VectorXd& whtscheb)
		{
			assert((n >= 1) && "n must be a positive integer");

			xscheb = ((Eigen::VectorXd::LinSpaced(n, 1.0, n).array() - 0.5) / n * M_PI).cos().reverse();
			whtscheb.setConstant(n, 2.0);

			double tj, tjm1, tjm2;
			for (Eigen::Index i = 0; i < n; ++i)
			{
				tjm2 = 1.0;
				tjm1 = xscheb(i);
				for (int j = 2; j < n; ++j)
				{
					tj = 2.0 * xscheb(i) * tjm1 - tjm2;
					tjm2 = tjm1;
					tjm1 = tj;
					if ((j % 2) == 0)
					{
						whtscheb(i) += 2.0 * tj * (1.0 / (j + 1.0) - 1.0 / (j - 1.0));
					}
				}
			}
			whtscheb *= 1.0 / n;
		}

		Eigen::VectorXd expansion_coefficients(
			const Eigen::Index n,
			const Eigen::VectorXd& vals)
		{
			assert((n >= 1) && "n must be a positive integer");

			Eigen::VectorXd xscheb, whtscheb;
			Eigen::MatrixXd umatr, vmatr;
			quadrature(n, xscheb, whtscheb);
			coefficient_matrix(n, xscheb, whtscheb, umatr, vmatr);
			return umatr * vals;
		}

		double evaluate_expansion(
			const Eigen::Index n,
			const Eigen::VectorXd& coefs,
			const double x)
		{
			assert((n >= 0) && "n must be a non-negative integer");

			if (n == 0)
			{
				return coefs(0);
			}
			double pj, pjm2, pjm1;
			pjm2 = 1.0;
			pjm1 = x;
			double val = coefs(0) * pjm2 + coefs(1) * pjm1;
			for (Eigen::Index j = 2; j <= n; ++j)
			{
				pj = 2.0 * x * pjm1 - pjm2;
				pjm2 = pjm1;
				pjm1 = pj;
				val += coefs(j) * pj;
			}
			return val;
		}

		void coefficient_matrix(
			const Eigen::Index n,
			const Eigen::VectorXd& xscheb,
			const Eigen::VectorXd& whtscheb,
			Eigen::MatrixXd& umatr,
			Eigen::MatrixXd& vmatr)
		{
			assert((n >= 1) && "n must be a positive integer");

			umatr.resize(n, n);
			vmatr.resize(n, n);

			double tj, tjm1, tjm2;
			for (Eigen::Index i = 0; i < n; ++i)
			{
				tjm2 = 1.0;
				tjm1 = xscheb(i);
				for (int j = 2; j < n; ++j)
				{
					tj = 2.0 * xscheb(i) * tjm1 - tjm2;
					tjm2 = tjm1;
					tjm1 = tj;
					vmatr(i, j) = tj;
				}
			}
			vmatr.col(0).setConstant(1.0 / n);
			vmatr.col(1) = xscheb.normalized() * std::sqrt(2.0 / n);
			vmatr.block(0, 2, n, n - 2).colwise().normalize();
			vmatr.block(0, 2, n, n - 2) *= std::sqrt(2.0 / n);
			umatr = vmatr.transpose();
			vmatr.col(0).setConstant(1.0);
			vmatr.block(0, 1, n, n - 1) *= 0.5 * n;
		}
	}

	namespace chebyshev_pract
	{
		void quadrature(
			const Eigen::Index n,
			Eigen::VectorXd& xscheb,
			Eigen::VectorXd& whtscheb)
		{
			assert((n >= 1) && "n must be a positive integer");

			xscheb = (Eigen::VectorXd::LinSpaced(n, 0.0, n - 1.0).array() / (n - 1.0) * M_PI).cos().reverse();
			whtscheb.setConstant(n, 1.0);

			double tj, tjm1, tjm2;
			for (Eigen::Index i = 0; i < n; ++i)
			{
				tjm2 = 1.0;
				tjm1 = xscheb(i);
				for (Eigen::Index j = 2; j < n; ++j)
				{
					tj = 2.0 * xscheb(i) * tjm1 - tjm2;
					tjm2 = tjm1;
					tjm1 = tj;
					if ((j % 2) == 0)
					{
						whtscheb(i) += tj * (1.0 / (j + 1.0) - 1.0 / (j - 1.0));
					}
				}
			}
			whtscheb *= 2.0 / (n - 1.0);
			whtscheb(0) *= 0.5;
			whtscheb(n - 1) *= 0.5;
		}

		void coefficient_matrix(
			const Eigen::Index n,
			const Eigen::VectorXd& xscheb,
			const Eigen::VectorXd& whtscheb,
			Eigen::MatrixXd& umatr,
			Eigen::MatrixXd& vmatr)
		{
			assert((n >= 1) && "n must be a positive integer");

			umatr.resize(n, n);
			vmatr.resize(n, n);

			double tj, tjm1, tjm2;
			for (Eigen::Index i = 0; i < n; ++i)
			{
				tjm2 = 1.0;
				tjm1 = xscheb(i);
				for (Eigen::Index j = 2; j < n; ++j)
				{
					tj = 2.0 * xscheb(i) * tjm1 - tjm2;
					tjm2 = tjm1;
					tjm1 = tj;
					vmatr(i, j) = tj;
				}
			}
			vmatr.col(0).setConstant(1.0);
			vmatr.col(1) = xscheb;
			umatr = vmatr.transpose();
			umatr *= 2.0 / (n - 1.0);
			umatr.col(0) *= 0.5;
			umatr.row(0) *= 0.5;
			umatr.col(n - 1) *= 0.5;
			umatr.row(n - 1) *= 0.5;
		}

		Eigen::MatrixXd differentiation_matrix(
			const Eigen::Index n,
			const Eigen::VectorXd& xscheb)
		{
			assert((n >= 1) && "n must be a positive integer");

			Eigen::MatrixXd dmatr(n, n);
			double ci = -1.0;
			double cj = -1.0;
			for (Eigen::Index i = 0; i < n; ++i)
			{
				ci = -ci;
				cj = -1.0;
				for (Eigen::Index j = 0; j < n; ++j)
				{
					cj = -cj;
					if (i == j)
					{
						dmatr(i, j) = -0.5 * xscheb(i) / (1.0 - std::pow(xscheb(i), 2.0));
					}
					else
					{
						dmatr(i, j) = (ci * cj) / (xscheb(i) - xscheb(j));
					}
				}
			}
			dmatr.col(0) *= 0.5;
			dmatr.row(0) *= 2.0;
			dmatr.col(n - 1) *= 0.5;
			dmatr.row(n - 1) *= 2.0;
			dmatr(0, 0) = -(1.0 / 6.0) * (1.0 + 2.0 * std::pow(n - 1.0, 2.0));
			dmatr(n - 1, n - 1) = -dmatr(0, 0);
			return dmatr;
		}
	}
}
