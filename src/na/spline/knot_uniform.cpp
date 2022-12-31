#include "na/spline/knot_uniform.h"
#include <cassert>
#include "Eigen/Core"
#include "na/core/vector.h"

namespace na
{
	namespace spline
	{
		Eigen::VectorXd open_knots(
			const Eigen::Index number,
			const Eigen::Index degree,
			const Eigen::Index regularity)
		{
			assert((regularity >= 0) && "open_knots: regularity must be a non-negative integer");
			assert((degree - regularity > 0) && "open_knots: the requested regularity is too high");

			Eigen::VectorXd knots(number * (degree - regularity) + 2 * regularity + 2);
			knots << Eigen::VectorXd::Zero(regularity + 1), na::repelem(Eigen::VectorXd::LinSpaced(number, 0.0, 1.0), degree - regularity), Eigen::VectorXd::Ones(regularity + 1);
			return knots;
		}

		Eigen::VectorXd periodic_knots(
			const Eigen::Index number,
			const Eigen::Index degree)
		{
			assert((number >= 2) && "periodic_knots: number must be at least 2");

			Eigen::VectorXd knots(number + 2 * degree);
			knots.segment(degree, number).setLinSpaced(0.0, 1.0);
			knots.head(degree) = (-1.0 / (number - 1.0)) * Eigen::VectorXi::LinSpaced(degree, 1, degree).cast<double>();
			knots.tail(degree) = (+1.0 / (number - 1.0)) * Eigen::VectorXi::LinSpaced(degree, 1, degree).cast<double>().array() + 1.0;
			return knots;
		}
	}
}
