#ifndef NA_SPLINE_NRB_REFINE_H
#define NA_SPLINE_NRB_REFINE_H
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity refine_curve(
			const NURBSEntity& nrb,
			const Eigen::Ref<const Eigen::ArrayXd>& x);
	}
}

#endif // !NA_SPLINE_NRB_REFINE_H
