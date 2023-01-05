#ifndef NA_INTEGRATE_DOPRI5_H
#define NA_INTEGRATE_DOPRI5_H
#include <functional>
#include "Eigen/Core"

namespace na
{
	void dopri5(
		const int n,
		const std::function<void(double, double*, double*)> f,
		double& x,
		const double xend,
		double* y,
		const double* rtol,
		const double* atol,
		const bool itol,
		int& status,
		int& nfcn,
		int& nstep,
		int& naccpt,
		int& nrejct,
		double safe = 0.0,
		double facl = 0.0,
		double facr = 0.0,
		double beta = 0.0,
		double hmax = 0.0,
		double h = 0.0,
		int nmax = 0,
		int nstiff = 0);
	
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
		const Eigen::VectorXd& y0,
		const double tolrel,
		const double tolabs,
		Eigen::VectorXd& y);

	bool dopri5(
		const std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> f,
		const double t0,
		const double t1,
		const Eigen::VectorXd& y0,
		const double tolrel,
		const double tolabs,
		Eigen::Ref<Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>> y);
}

#endif // !NA_INTEGRATE_DOPRI5_H
