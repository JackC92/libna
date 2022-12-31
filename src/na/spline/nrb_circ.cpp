#include "na/spline/nrb_circ.h"
#include <cmath>
#include "Eigen/Core"
#include "na/geometry/transformation.h"
#include "na/spline/knot_uniform.h"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity nrb_circ(
			const double radius,
			const Eigen::Vector3d& center,
			const double start_angle,
			const double end_angle)
		{
			double sweep = end_angle - start_angle;
			if (sweep < 0.0)
			{
				sweep += 2.0 * M_PI;
			}

			Eigen::Index narcs;
			Eigen::VectorXd knots;
			if (std::abs(sweep) <= M_PI_2)
			{
				narcs = 1;
				knots = na::spline::open_knots(2, 2, 0);
			}
			else if (std::abs(sweep) <= M_PI)
			{
				narcs = 2;
				knots = na::spline::open_knots(3, 2, 0);
			}
			else if (std::abs(sweep) <= 1.5 * M_PI)
			{
				narcs = 3;
				knots = na::spline::open_knots(4, 2, 0);
			}
			else
			{
				narcs = 4;
				knots = na::spline::open_knots(5, 2, 0);
			}

			double dsweep = sweep / (2.0 * narcs);
			double wm = std::cos(dsweep);
			double x = radius * wm;
			double y = radius * std::sin(dsweep);
			double xm = x + y * std::tan(dsweep);

			Eigen::MatrixX4d arc_coefs(3, 4), coefs(2 * narcs + 1, 4);
			arc_coefs << x, -y, 0.0, 1.0,
						 wm * xm, 0.0, 0.0, wm,
						 x, y, 0.0, 1.0;
			coefs.block<3, 4>(0, 0) = arc_coefs * na::homogeneous_rotation_angleZ(start_angle + dsweep).transpose();
			Eigen::Matrix4d R = na::homogeneous_rotation_angleZ(2.0 * dsweep);
			for (Eigen::Index n = 1; n < narcs; ++n)
			{
				coefs.middleRows(2 * n + 1, 2) = coefs.middleRows(2 * n - 1, 2) * R.transpose();
			}
			coefs = coefs * na::homogeneous_translation(center).transpose();
			return NURBSEntity(2, 2 * narcs + 1, knots, coefs, false);
		}
	}
}
