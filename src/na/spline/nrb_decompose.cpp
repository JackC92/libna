#include "na/spline/nrb_decompose.h"
#include <vector>
#include "Eigen/Core"
#include "na/core/array.h"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		std::vector<Eigen::MatrixX4d> decompose_curve(const NURBSEntity& nrb)
		{
			Eigen::Index deg = nrb.degree()[0];
			Eigen::Index size;
			{
				Eigen::ArrayXd knots = nrb.knots()[0];
				na::deduplicate<double>(knots);
				size = knots.size();
			}
			std::vector<Eigen::MatrixX4d> curves;
			curves.reserve(size - 1);

			Eigen::Index m = nrb.number()[0] + deg;
			Eigen::Index a = deg;
			Eigen::Index b = deg + 1;
			Eigen::Index nb = 0;
			curves.emplace_back(nrb.coefs().topRows(deg + 1));

			Eigen::Ref<const Eigen::ArrayXd> knots(nrb.knots()[0]);
			Eigen::ArrayXd alphas(deg - 1);
			while (b < m)
			{
				Eigen::Index i = b;
				while ((b < m) && (knots.coeffRef(b + 1) == knots.coeffRef(b)))
				{
					b += 1;
				}
				if (b < m)
				{
					curves.emplace_back(Eigen::MatrixX4d::Zero(deg + 1, 4));
				}
				Eigen::Index mult = b - i + 1;
				if (mult < deg)
				{
					double numer = knots.coeffRef(b) - knots.coeffRef(a);
					for (Eigen::Index j = deg; j > mult; --j)
					{
						alphas.coeffRef(j - mult - 1) = numer / (knots.coeffRef(a + j) - knots.coeffRef(a));
					}
					Eigen::Index r = deg - mult;
					for (Eigen::Index j = 1; j <= r; ++j)
					{
						Eigen::Index save = r - j;
						Eigen::Index s = deg + j;
						for (Eigen::Index k = deg; k >= s; --k)
						{
							double alpha = alphas.coeffRef(k - s);
							curves[nb].row(k) = alpha * curves[nb].row(k) + (1.0 - alpha) * curves[nb].row(k - 1);
						}
						if (b < m)
						{
							curves[nb + 1].row(save) = curves[nb].row(deg);
						}
					}
				}
				nb += 1;
				if (b < m)
				{
					for (Eigen::Index i = deg - mult; i <= deg; ++i)
					{
						curves[nb].row(i) = nrb.coefs().row(b - deg + i);
					}
					a = b;
					b = b + 1;
				}
			}
			return curves;
		}
	}
}