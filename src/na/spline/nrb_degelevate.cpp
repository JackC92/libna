#include "na/spline/nrb_degelevate.h"
#include "Eigen/Core"
#include "na/spline/bsp_degelevate.h"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity nrb_degelevate_curve(
			const NURBSEntity& nrb, 
			const Eigen::Index t)
		{
			Eigen::ArrayXd knots;
			Eigen::MatrixXd coefs;
			na::bspline::bsp_degelevate(nrb.coefs(), nrb.knots()[0], nrb.degree()[0], t, coefs, knots);
			return NURBSEntity(nrb.degree()[0] + t, coefs.rows(), knots, coefs);
		}
	}
}
