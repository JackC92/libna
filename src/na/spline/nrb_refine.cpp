#include "na/spline/nrb_refine.h"
#include "Eigen/Core"
#include "na/spline/bsp_refine.h"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity refine_curve(const NURBSEntity& nrb, const Eigen::Ref<const Eigen::ArrayXd>& x)
		{
			Eigen::ArrayXd knots;
			Eigen::MatrixXd coefs;
			na::bspline::refine_curve(nrb.coefs(), nrb.knots()[0], nrb.degree()[0], x, coefs, knots);
			return NURBSEntity(nrb.degree()[0], coefs.rows(), knots, coefs);
		}
	}
}
