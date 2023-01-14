#ifndef NA_SPLINE_FIND_SPAN_H
#define NA_SPLINE_FIND_SPAN_H
#include "Eigen/Core"

namespace na
{
	namespace spline
	{
		Eigen::Index find_span(
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index degree,
			const double u);

		Eigen::ArrayXi find_span(
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index degree,
			const Eigen::Ref<const Eigen::ArrayXd>& u);
	}
}

#endif // !NA_SPLINE_FIND_SPAN_H
