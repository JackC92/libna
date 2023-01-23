#include "na/spline/bsp_refine.h"
#include <cmath>
#include "Eigen/Core"
#include "na/spline/find_span.h"

namespace na
{
	namespace bspline
	{
		void refine_curve(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index deg,
			const Eigen::Ref<const Eigen::ArrayXd>& x,
			Eigen::MatrixXd& coefsout,
			Eigen::ArrayXd& knotsout)
		{
			if (x.size() == 0)
			{
				coefsout = coefs;
				knotsout = knots;
				return;
			}

			knotsout.resize(knots.size() + x.size());
			coefsout.resize(coefs.rows() + x.size(), coefs.cols());
			Eigen::Index m = knots.size() - 1;
			Eigen::Index n = coefs.rows() - 1;
			Eigen::Index r = x.size() - 1;
			Eigen::Index a = na::spline::find_span(knots, deg, x.coeffRef(0));
			Eigen::Index b = na::spline::find_span(knots, deg, x.coeffRef(x.size() - 1));
			b += 1;
			coefsout.topRows(a - deg) = coefs.topRows(a - deg);
			coefsout.bottomRows(n - b + 2) = coefs.bottomRows(n - b + 2);
			knotsout.head(a + 1) = knots.head(a + 1);
			knotsout.tail(m - b - deg + 1) = knots.tail(m - b - deg + 1);
			Eigen::Index i = b + deg - 1;
			Eigen::Index k = b + deg + r;
			for (Eigen::Index j = r; j >= 0; --j)
			{
				while ((x.coeffRef(j) <= knots.coeffRef(i)) && (i > a))
				{
					coefsout.row(k - deg - 1) = coefs.row(i - deg - 1);
					knotsout.coeffRef(k) = knots.coeffRef(i);
					k -= 1;
					i -= 1;
				}
				coefsout.row(k - deg - 1) = coefsout.row(k - deg);
				for (Eigen::Index l = 1; l <= deg; ++l)
				{
					Eigen::Index ind = k - deg + l;
					double alfa = knotsout.coeffRef(k + l) - x.coeffRef(j);
					if (std::abs(alfa) == 0.0)
					{
						coefsout.row(ind - 1) = coefsout.row(ind);
					}
					else
					{
						alfa /= knotsout.coeffRef(k + l) - knots.coeffRef(i - deg + l);
						coefsout.row(ind - 1) = alfa * coefsout.row(ind - 1) + (1.0 - alfa) * coefsout.row(ind);
					}
				}
				knotsout.coeffRef(k) = x.coeffRef(j);
				k -= 1;
			}
		}
	}
}
