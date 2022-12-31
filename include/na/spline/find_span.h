#ifndef NA_SPLINE_FIND_SPAN_H
#define NA_SPLINE_FIND_SPAN_H
#include <vector>
#include "Eigen/Core"

namespace na
{
	namespace spline
	{
		Eigen::Index find_span(
			Eigen::Ref<const Eigen::VectorXd> knots,
			const Eigen::Index degree,
			const double u);

		std::vector<Eigen::Index> find_span(
			Eigen::Ref<const Eigen::VectorXd> knots,
			const Eigen::Index degree,
			Eigen::Ref<const Eigen::VectorXd> u);
	}
}

#endif // !NA_SPLINE_FIND_SPAN_H
