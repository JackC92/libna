#ifndef NA_SPLINE_NRB_UNCLAMP_H
#define NA_SPLINE_NRB_UNCLAMP_H
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity unclamp_curve(
			const NURBSEntity& nrb,
			const Eigen::Index k = -1);
	}
}

#endif // !NA_SPLINE_NRB_UNCLAMP_H
