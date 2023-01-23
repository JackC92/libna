#include "na/spline/nrb_unclamp.h"
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity unclamp_curve(
			const NURBSEntity& nrb,
			const Eigen::Index k)
		{
			Eigen::Index deg = nrb.degree()[0];
			Eigen::Index n = nrb.coefs().rows();
			Eigen::Index m = nrb.knots()[0].size();
			Eigen::Index kk = k;
			if ((kk < 0) || (kk >= deg))
			{
				kk = deg - 1;
			}
			Eigen::ArrayXd knots = nrb.knots()[0];
			Eigen::MatrixXd coefs = nrb.coefs();
			for (Eigen::Index ii = 0; ii <= kk; ++ii)
			{
				knots.coeffRef(kk - ii) = knots.coeffRef(kk - ii + 1) - (knots.coeffRef(n - ii) - knots.coeffRef(n - ii - 1));
			}
			for (Eigen::Index ii = deg - kk - 1; ii <= deg - 2; ++ii)
			{
				for (Eigen::Index jj = ii; jj >= 0; --jj)
				{
					double alpha = (knots.coeffRef(deg + jj + 1) - knots.coeffRef(deg + jj - ii - 1)) / (knots.coeffRef(deg + jj + 1) - knots.coeffRef(deg));
					coefs.row(jj) = alpha * coefs.row(jj) + (1.0 - alpha) * coefs.row(jj + 1);
				}
			}
			for (Eigen::Index ii = 0; ii <= kk; ++ii)
			{
				knots.coeffRef(m - kk + ii - 1) = knots.coeffRef(m - kk + ii - 2) + knots.coeffRef(deg + ii + 1) - knots.coeffRef(deg + ii);
			}
			for (Eigen::Index ii = deg - kk; ii <= deg - 1; ++ii)
			{
				for (Eigen::Index jj = ii; jj > 0; --jj)
				{
					double alpha = (knots.coeffRef(n - jj + ii + 2) - knots.coeffRef(n - jj)) / (knots.coeffRef(n) - knots.coeffRef(n - jj));
					coefs.row(n - jj) = alpha * coefs.row(n - jj) + (1.0 - alpha) * coefs.row(n - jj - 1);
				}
			}
			return NURBSEntity(deg, n, knots, coefs);
		}
	}
}
