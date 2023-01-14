#ifndef NA_SPLINE_BSP_DEGELEVATE_H
#define NA_SPLINE_BSP_DEGELEVATE_H
#include "Eigen/Core"

namespace na
{
	namespace bspline
	{
		void bsp_degelevate(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index degree,
			const Eigen::Index t,
			Eigen::MatrixXd& coefsout,
			Eigen::ArrayXd& knotsout);
	}
}

#endif // !NA_SPLINE_BSP_DEGELEVATE_H
