#ifndef NA_SPLINE_BSP_KNTREFINE_H
#define NA_SPLINE_BSP_KNTREFINE_H
#include "Eigen/Core"

namespace na
{
	namespace bspline
	{
		void kntrefine_curve(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index degree,
			const Eigen::Ref<const Eigen::ArrayXd>& x,
			Eigen::MatrixXd& coefsout,
			Eigen::ArrayXd& knotsout);
	}
}

#endif // !NA_SPLINE_BSP_KNTREFINE_H
