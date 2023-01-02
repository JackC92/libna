#ifndef NA_INTEGRATE_RUNGE_KUTTA_H
#define NA_INTEGRATE_RUNGE_KUTTA_H
#include "Eigen/Core"

namespace na
{
	bool dopri5(
		const std::function<double(double, double)> f,
		const double t0,
		const double t1,
		const double y0,
		const double tolrel,
		const double tolabs,
		double& y);
	
	bool dopri5(
		const std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> f,
		const double t0,
		const double t1,
		const Eigen::Ref<const Eigen::VectorXd>& y0,
		const double tolrel,
		const double tolabs,
		Eigen::Ref<Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>> y);
}

#endif // !NA_INTEGRATE_RUNGE_KUTTA_H
