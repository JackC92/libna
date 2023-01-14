#ifndef NA_SPLINE_NRB_CIRC_H
#define NA_SPLINE_NRB_CIRC_H
#include <cmath>
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity nrb_circ(
			const double radius = 1.0,
			const Eigen::Vector3d& center = Eigen::Vector3d::Zero(),
			const double start_angle = 0.0,
			const double end_angle = 2.0 * M_PI);
	}
}

#endif // !NA_SPLINE_NRB_CIRC_H
