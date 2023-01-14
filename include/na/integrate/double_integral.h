#ifndef NA_INTEGRATE_DOUBLE_INTEGRAL_H
#define NA_INTEGRATE_DOUBLE_INTEGRAL_H
#include <cmath>
#include <functional>
#include "Eigen/Core"
#include "na/special/legendre.h"
#include "na/type_traits/is_array_or_matrix.h"
#include "na/type_traits/is_float_or_complex.h"

namespace na
{
	// Compute the integral of f(x, y) over the disk of radius R centered at the origin
	template <
		typename ReturnType,
		std::enable_if_t<na::is_float_or_complex_v<ReturnType>, bool> = true>
	ReturnType integrate_disk(
		const std::function<ReturnType(double, double)> f,
		const double R)
	{
		Eigen::VectorXd xslege, whtslege;
		na::legendre::quadrature(4, xslege, whtslege);

		double val = 0.0;
		// Compute the integral of the function in each quarter of the domain 0 <= theta < 2 * pi.
		// This function is included to have comparable results with section 6.1 of the paper
		//   Myung-Jin Choi, Roger A. Sauer, and Sven Klinkel, "An isogeometric finite element formulation for geometrically exact Timoshenko beams with extensible directors",
		//   Computer Methods in Applied Mechanics and Engineering, 2021.
		Eigen::MatrixXd fvals(xslege.size(), xslege.size());
		for (int i = 0; i < 4; ++i)
		{
			double theta1 = M_PI_2 * i;
			double theta2 = M_PI_2 * (i + 1.0);
			for (Eigen::Index ir = 0; ir < xslege.size(); ++ir)
			{
				// Apply the Gauss-Legendre quadrature rule to the inner integral after a change of interval to [0, 1].
				double r = 0.5 * (xslege(ir) + 1.0);
				for (Eigen::Index it = 0; it < xslege.size(); ++it)
				{
					double theta = 0.5 * ((theta2 - theta1) * xslege(it) + (theta2 + theta1));

					double x = R * r * std::cos(theta);
					double y = R * r * std::sin(theta);
					fvals(it, ir) = f(x, y) * r;
				}
			}
			val += (theta2 - theta1) * whtslege.transpose() * fvals * whtslege;
		}
		return 0.25 * R * R * val;
	}

	template <
		typename ReturnType,
		std::enable_if_t<na::is_array_or_matrix_v<ReturnType>, bool> = true>
	ReturnType integrate_disk(
		const std::function<ReturnType(double, double)> f,
		const double R)
	{
		Eigen::VectorXd xslege, whtslege;
		na::legendre::quadrature(4, xslege, whtslege);

		ReturnType val = ReturnType::Zero(f(0.0, 0.0).rows(), f(0.0, 0.0).cols());
		for (int i = 0; i < 4; ++i)
		{
			double theta1 = M_PI_2 * i;
			double theta2 = M_PI_2 * (i + 1.0);
			for (Eigen::Index ir = 0; ir < xslege.size(); ++ir)
			{
				// Apply the Gauss-Legendre quadrature rule to the inner integral after a change of interval to [0, 1].
				double r = 0.5 * (xslege(ir) + 1.0);
				for (Eigen::Index it = 0; it < xslege.size(); ++it)
				{
					double theta = 0.5 * ((theta2 - theta1) * xslege(it) + (theta2 + theta1));

					double x = R * r * std::cos(theta);
					double y = R * r * std::sin(theta);
					val += (theta2 - theta1) * whtslege(ir) * whtslege(it) * f(x, y) * r;
				}
			}
		}
		return 0.25 * R * R * val;
	}

	// Compute the integral of f(x, y) over the rectangle with bottom-left corner (x1, y1) and top-right corner (x2, y2)
	template <
		typename ReturnType,
		std::enable_if_t<na::is_float_or_complex_v<ReturnType>, bool> = true>
	ReturnType integrate_rectangle(
		const std::function<double(double, double)> f,
		const double x1,
		const double x2,
		const double y1,
		const double y2)
	{
		Eigen::VectorXd xslege, whtslege;
		na::legendre::quadrature(4, xslege, whtslege);

		double val = 0.0;
		Eigen::MatrixXd fvals(xslege.size(), xslege.size());
		for (Eigen::Index iy = 0; iy < xslege.size(); ++iy)
		{
			for (Eigen::Index ix = 0; ix < xslege.size(); ++ix)
			{
				double x = 0.5 * ((x2 - x1) * xslege(ix) + (x2 + x1));
				double y = 0.5 * ((y2 - y1) * xslege(iy) + (y2 + y1));
				fvals(ix, iy) = f(x, y);
			}
		}
		return 0.25 * (x2 - x1) * (y2 - y1) * whtslege.transpose() * fvals * whtslege;
	}

	template <
		typename ReturnType,
		std::enable_if_t<na::is_array_or_matrix_v<ReturnType>, bool> = true>
	ReturnType integrate_rectangle(
		const std::function<ReturnType(double, double)> f,
		const double x1,
		const double x2,
		const double y1,
		const double y2)
	{
		Eigen::VectorXd xslege, whtslege;
		na::legendre::quadrature(4, xslege, whtslege);

		ReturnType val = ReturnType::Zero(f(0.0, 0.0).rows(), f(0.0, 0.0).cols());
		for (Eigen::Index ix = 0; ix < xslege.size(); ++ix)
		{
			for (Eigen::Index iy = 0; iy < xslege.size(); ++iy)
			{
				double x = 0.5 * ((x2 - x1) * xslege(ix) + (x2 + x1));
				double y = 0.5 * ((y2 - y1) * xslege(iy) + (y2 + y1));
				val += whtslege(ix) * whtslege(iy) * f(x, y);
			}
		}
		return 0.25 * (x2 - x1) * (y2 - y1) * val;
	}
}

#endif // !NA_INTEGRATE_DOUBLE_INTEGRAL_H
