#ifndef NA_SPLINE_NRB_CIRC_H
#define NA_SPLINE_NRB_CIRC_H
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity nrb_circ(
			const double radius,
			const Eigen::Vector3d& center,
			const double start_angle,
			const double end_angle);
	}
}

#endif // !NA_SPLINE_NRB_CIRC_H
