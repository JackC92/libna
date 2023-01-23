#ifndef NA_SPLINE_NRB_TESSELLATE_H
#define NA_SPLINE_NRB_TESSELLATE_H
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		void tessellate(
			const NURBSEntity& nrb,
			Eigen::MatrixX3d& V,
			const Eigen::Index upts = 50);

		void tessellate(
			const NURBSEntity& nrb,
			Eigen::MatrixX3d& V,
			Eigen::MatrixX3i& F,
			const Eigen::Index upts = 50,
			const Eigen::Index vpts = 50);
	}
}

#endif // !NA_SPLINE_NRB_TESSELLATE_H
