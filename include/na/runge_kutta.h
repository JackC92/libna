#ifndef NA_RUNGE_KUTTA_H
#define NA_RUNGE_KUTTA_H
#include <functional>

namespace na
{
	void dopri5(
		const std::function<double(double, double)> f,
		const double x0,
		const double x1,
		const double y0,
		const double tolrel,
		const double tolabs,
		double& x,
		double& y,
		double& h,
		int& state);
}

#endif // !NA_RUNGE_KUTTA_H
