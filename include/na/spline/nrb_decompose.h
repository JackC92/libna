#ifndef NA_SPLINE_NRB_DECOMPOSE_H
#define NA_SPLINE_NRB_DECOMPOSE_H
#include <vector>
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		// Decompose the input NURBS curve into Bezier segments and return only the control points of the Bezier segments.
		// For the returned control points of Bezier segments, the degree of the Bezier segment can be inferred from the number of rows.
		std::vector<Eigen::MatrixX4d> nrb_decompose_curve(const NURBSEntity& nrb);
	}
}

#endif // !NA_SPLINE_NRB_DECOMPOSE_H
