#ifndef NA_SPLINE_NRB_KNTREMOVE_H
#define NA_SPLINE_NRB_KNTREMOVE_H
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		Eigen::Index kntremove_curve(
			NURBSEntity& nrb,
			const Eigen::Index r,
			const Eigen::Index num,
			const double tol = 0.0);
	}
}

#endif // !NA_SPLINE_NRB_KNTREMOVE_H
