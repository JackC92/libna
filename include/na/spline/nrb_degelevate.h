#ifndef NA_SPLINE_NRB_DEGELEVATE_H
#define NA_SPLINE_NRB_DEGELEVATE_H
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity degelevate_curve(
			const NURBSEntity& nrb,
			const Eigen::Index t);
	}
}

#endif // !NA_SPLINE_NRB_DEGELEVATE_H
