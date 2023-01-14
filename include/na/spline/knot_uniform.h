#ifndef NA_SPLINE_KNOT_UNIFORM_H
#define NA_SPLINE_KNOT_UNIFORM_H
#include "Eigen/Core"

namespace na
{
	namespace spline
	{
		// Generate an open uniform knot vector spanning the interval [0, 1]
		//
		// Input parameters:
		//   number: number of distinct knots in [0, 1]
		//   degree: degree of the knot vector
		//   regularity: the regularity at interior knots
		//
		// Output parameters:
		//   knots: an open uniform knot vector
		Eigen::ArrayXd open_knots(
			const Eigen::Index number,
			const Eigen::Index degree,
			const Eigen::Index regularity);

		// Generate a periodic uniform knot vector spanning the interval [0, 1]
		//
		// Input parameters:
		//   number: number of distinct knots in [0, 1]
		//   degree: degree of the knot vector
		//
		// Output parameters:
		//   knots: a periodic uniform knot vector
		Eigen::ArrayXd periodic_knots(
			const Eigen::Index number,
			const Eigen::Index degree);
	}
}

#endif // !NA_SPLINE_KNOT_UNIFORM_H
