#include "na/spline/find_span.h"
#include <cassert>
#include "Eigen/Core"

namespace na
{
	namespace spline
	{
		Eigen::Index find_span(
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index deg,
			const double u)
		{
			assert(((u >= knots.coeffRef(deg)) && (u <= *(knots.end() - deg - 1))) && "find_span: u is out of range");

			Eigen::Index n = knots.size() - deg - 2;
			if (u == knots.coeffRef(n + 1))
			{
				return n;
			}
			Eigen::Index low = deg;
			Eigen::Index high = n + 1;
			Eigen::Index mid = (low + high) / 2;
			while (u < knots.coeffRef(mid) || u >= knots.coeffRef(mid + 1))
			{
				if (u < knots.coeffRef(mid))
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
		
		Eigen::ArrayXi find_span(
			const Eigen::Ref<const Eigen::ArrayXd>& knots, 
			const Eigen::Index deg,
			const Eigen::Ref<const Eigen::ArrayXd>& u)
		{
			Eigen::ArrayXi span(u.size());
			for (Eigen::Index i = 0; i < u.size(); ++i)
			{
				span.coeffRef(i) = find_span(knots, deg, u.coeffRef(i));
			}
			return span;
		}
	}
}
