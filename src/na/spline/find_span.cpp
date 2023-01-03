#include "na/spline/find_span.h"
#include <vector>
#include "Eigen/Core"

namespace na
{
	namespace spline
	{
		Eigen::Index find_span(
			const Eigen::Ref<const Eigen::VectorXd>& knots,
			const Eigen::Index degree,
			const double u)
		{
			assert(((u >= knots(degree)) && (u <= *(knots.end() - degree - 1))) && "find_span: u is out of range");

			Eigen::Index n = knots.size() - degree - 1;
			if (u == knots(n))
			{
				return n - 1;
			}
			Eigen::Index low = degree;
			Eigen::Index high = n;
			Eigen::Index mid = (low + high) / 2;
			while (u < knots(mid) || u >= knots(mid + 1))
			{
				if (u < knots(mid))
				{
					high = mid;
				}
				else
				{
					low = mid;
				}
				mid = (low + high) / 2;
			}
			return mid;
		}

		std::vector<Eigen::Index> find_span(
			const Eigen::Ref<const Eigen::VectorXd>& knots,
			const Eigen::Index degree,
			const Eigen::Ref<const Eigen::VectorXd>& u)
		{
			std::vector<Eigen::Index> indices;
			indices.reserve(u.size());
			for (Eigen::Index i = 0; i < u.size(); ++i)
			{
				indices.emplace(indices.begin() + i, find_span(knots, degree, u(i)));
			}
			return indices;
		}
	}
}
