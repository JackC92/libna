#include "na/spline/knot_uniform.h"
#include <cassert>
#include "Eigen/Core"
#include "na/core/array.h"

namespace na
{
	namespace spline
	{
		Eigen::ArrayXd open_knots(
			const Eigen::Index number,
			const Eigen::Index degree,
			const Eigen::Index regularity)
		{
			assert((regularity >= 0) && "open_knots: regularity must be a non-negative integer");
			assert((degree - regularity > 0) && "open_knots: the requested regularity is too high");

			Eigen::ArrayXd knots(number * (degree - regularity) + 2 * regularity + 2);
			knots << Eigen::ArrayXd::Zero(regularity + 1), na::repelem(Eigen::ArrayXd::LinSpaced(number, 0.0, 1.0), degree - regularity), Eigen::ArrayXd::Ones(regularity + 1);
			return knots;
		}

		Eigen::ArrayXd periodic_knots(
			const Eigen::Index number,
			const Eigen::Index degree)
		{
			assert((number >= 2) && "periodic_knots: number must be at least 2");

			Eigen::ArrayXd knots(number + 2 * degree);
			knots.segment(degree, number).setLinSpaced(0.0, 1.0);
			knots.head(degree) = (-1.0 / (number - 1.0)) * Eigen::ArrayXd::LinSpaced(degree, 1, degree);
			knots.tail(degree) = (+1.0 / (number - 1.0)) * Eigen::ArrayXd::LinSpaced(degree, 1, degree) + 1.0;
			return knots;
		}
	}
}
