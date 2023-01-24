#include "na/spline/nrb_kntrefine.h"
#include "Eigen/Core"
#include "na/spline/bsp_kntrefine.h"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity kntrefine_curve(const NURBSEntity& nrb, const Eigen::Ref<const Eigen::ArrayXd>& x)
		{
			Eigen::ArrayXd knots;
			Eigen::MatrixXd coefs;
			na::bspline::kntrefine_curve(nrb.coefs(), nrb.knots()[0], nrb.degree()[0], x, coefs, knots);
			return NURBSEntity(nrb.degree()[0], coefs.rows(), knots, coefs);
		}
	}
}
