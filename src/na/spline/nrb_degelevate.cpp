#include "na/spline/nrb_degelevate.h"
#include "Eigen/Core"
#include "na/spline/bsp_degelevate.h"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity degelevate_curve(
			const NURBSEntity& nrb, 
			const Eigen::Index t)
		{
			Eigen::ArrayXd knots;
			Eigen::MatrixXd coefs;
			na::bspline::degelevate(nrb.coefs(), nrb.knots()[0], nrb.degree()[0], t, coefs, knots);
			return NURBSEntity(nrb.degree()[0] + t, coefs.rows(), knots, coefs);
		}
	}
}
