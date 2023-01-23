#ifndef NA_SPLINE_BSP_REFINE_H
#define NA_SPLINE_BSP_REFINE_H
#include "Eigen/Core"

namespace na
{
	namespace bspline
	{
		void refine_curve(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index degree,
			const Eigen::Ref<const Eigen::ArrayXd>& x,
			Eigen::MatrixXd& coefsout,
			Eigen::ArrayXd& knotsout);
	}
}

#endif // !NA_SPLINE_BSP_REFINE_H
