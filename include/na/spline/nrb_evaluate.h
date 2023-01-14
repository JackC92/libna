#ifndef NA_SPLINE_NRB_EVALUATE_H
#define NA_SPLINE_NRB_EVALUATE_H
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		// Evaluate the NURBS curve at the given parametric coordinate u.
		Eigen::Vector3d evaluate(
			const NURBSEntity& nrb,
			const double u);

		// Evaluate the derivatives of the NURBS curve at the given parametric coordinate u.
		Eigen::MatrixX3d evaluate(
			const NURBSEntity& nrb,
			const Eigen::Index m,
			const double u);

		// Evaluate the NURBS surface at the given parametric coordinate (u, v).
		Eigen::Vector3d evaluate(
			const NURBSEntity& nrb,
			const double u,
			const double v);

		// Evaluate the derivatives of the NURBS surface at the given parametric coordinate (u, v).
		// Return the derivatives of all order up to m, so there are (m + 1) * (m + 2) / 2 derivatives.
		Eigen::MatrixX3d evaluate(
			const NURBSEntity& nrb,
			const Eigen::Index m,
			const double u,
			const double v);
	}
}

#endif // !NA_SPLINE_NRB_EVALUATE_H
