#include "na/spline/nrb_line.h"
#include "Eigen/Core"
#include "na/spline/knot_uniform.h"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity nrb_line(
			const Eigen::Vector2d& start,
			const Eigen::Vector2d& end)
		{
			Eigen::VectorXd knots = na::spline::open_knots(2, 1, 0);
			Eigen::MatrixXd coefs(2, 2);
			coefs.row(0) = start;
			coefs.row(1) = end;
			return NURBSEntity(1, 2, knots, coefs, false);
		}

		NURBSEntity nrb_line(
			const Eigen::Vector3d& start,
			const Eigen::Vector3d& end)
		{
			Eigen::VectorXd knots = na::spline::open_knots(2, 1, 0);
			Eigen::MatrixXd coefs(2, 3);
			coefs.row(0) = start;
			coefs.row(1) = end;
			return NURBSEntity(1, 2, knots, coefs, false);
		}
	}
}
