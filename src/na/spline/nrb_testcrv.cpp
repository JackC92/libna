#include "na/spline/nrb_testcrv.h"
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity nrb_testcrv()
		{
			Eigen::ArrayXd knots(10);
			Eigen::MatrixX3d coefs(7, 3);
			knots << 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0;
			coefs << 0.5, 3.0, 0.0,
				1.5, 5.5, 0.0,
				4.5, 5.5, 0.0,
				3.0, 1.5, 0.0,
				7.5, 1.5, 0.0,
				6.0, 4.0, 0.0,
				8.5, 4.5, 0.0;
			return NURBSEntity(2, 7, knots, coefs);
		}
	}
}
