#ifndef NA_SPLINE_NRB_LINE_H
#define NA_SPLINE_NRB_LINE_H
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity nrb_line(
			const Eigen::Vector2d& start,
			const Eigen::Vector2d& end);

		NURBSEntity nrb_line(
			const Eigen::Vector3d& start,
			const Eigen::Vector3d& end);
	}
}

#endif // !NA_SPLINE_NRB_LINE_H
