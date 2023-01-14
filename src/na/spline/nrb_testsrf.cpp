#include "na/spline/nrb_testsrf.h"
#include "Eigen/Core"
#include "na/spline/knot_uniform.h"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		NURBSEntity nrb_testsrf()
		{
			Eigen::MatrixX3d coefs(25, 3);
			coefs << 0.0, 0.0, 2.0,
				3.0, 0.0, 2.0,
				5.0, 0.0, 7.0,
				8.0, 0.0, 7.0,
				10.0, 0.0, 8.0,
				0.0, 3.0, 0.0,
				3.0, 3.0, 0.0,
				5.0, 3.0, 5.0,
				8.0, 3.0, 5.0,
				10.0, 3.0, 7.0,
				0.0, 5.0, 0.0,
				3.0, 5.0, 0.0,
				5.0, 5.0, 5.0,
				8.0, 5.0, 5.0,
				10.0, 5.0, 7.0,
				0.0, 8.0, 5.0,
				3.0, 8.0, 5.0,
				5.0, 8.0, 8.0,
				8.0, 8.0, 8.0,
				10.0, 8.0, 10.0,
				0.0, 10.0, 5.0,
				3.0, 10.0, 5.0,
				5.0, 10.0, 8.0,
				8.0, 10.0, 8.0,
				10.0, 10.0, 10.0;
			return NURBSEntity({ 2, 2 }, { 5, 5 }, { na::spline::open_knots(4, 2, 1), na::spline::open_knots(4, 2, 1) }, coefs);
		}
	}
}
