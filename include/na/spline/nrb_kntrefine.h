#ifndef NA_SPLINE_NRB_KNTREFINE_H
#define NA_SPLINE_NRB_KNTREFINE_H
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity kntrefine_curve(
			const NURBSEntity& nrb,
			const Eigen::Ref<const Eigen::ArrayXd>& x);
	}
}

#endif // !NA_SPLINE_NRB_KNTREFINE_H
