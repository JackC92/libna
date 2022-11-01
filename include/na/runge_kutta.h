#ifndef NA_RUNGE_KUTTA_H
#define NA_RUNGE_KUTTA_H
#include <functional>
#include "Eigen/Core"

namespace na
{
	namespace internal
	{
		void dopri5_hinit(
			const std::function<double(double, double)> f,
			const double x0,
			const double x1,
			const double y0,
			const double tolrel,
			const double tolabs,
			double& h);

		template <typename DerivedIn, int Size, int Options = 0>
		void dopri5_hinit(
			const std::function<Eigen::Matrix<double, Size, 1, Options, Size, 1>(double, const Eigen::Matrix<double, Size, 1, Options, Size, 1>&)> f,
			const double x0,
			const double x1,
			const Eigen::MatrixBase<DerivedIn>& y0,
			const double tolrel,
			const double tolabs,
			double& h)
		{
			typedef typename Eigen::Matrix<double, Size, 1, Options, Size, 1> Vector;
			// Compute a first guess for explicit Euler as
			//   h = 0.01 * norm(y0) / norm(f(x0, y0))
			Vector f0 = f(x0, y0);
			Vector sk = tolabs + tolrel * y0.array().abs();
			double dnf = f0.cwiseQuotient(sk).squaredNorm();
			double dny = y0.cwiseQuotient(sk).squaredNorm();
			if ((dnf < 1e-10) || (dny < 1e-10))
			{
				h = 1e-6;
			}
			else
			{
				h = 0.01 * std::sqrt(dny / dnf);
			}
			// TODO: allow user-specified maximal step size
			h = std::min<double>(h, x1 - x0);
			Vector f1 = f(x0 + h, y0 + h * f0);
			double der2 = (f1 - f0).cwiseQuotient(sk).norm() / h;
			// Step size is computed such that
			//   h**5 * max(norm(f(x0, y0)), norm(der2)) = 0.01
			double der12 = std::max<double>(der2, std::sqrt(dnf));
			if (der12 < 1e-15)
			{
				h = std::min<double>({ 100.0 * h, std::max<double>(1e-6, h * 1e-3), x1 - x0 });
			}
			else
			{
				h = std::min<double>({ 100.0 * h, std::pow(0.01 / der12, 0.2), x1 - x0 });
			}
		}
	}

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

	template <typename DerivedIn, typename DerivedOut, int Size, int Options = 0>
	void dopri5(
		const std::function<Eigen::Matrix<double, Size, 1, Options, Size, 1>(double, const Eigen::Matrix<double, Size, 1, Options, Size, 1>&)> f,
		const double x0,
		const double x1,
		const Eigen::MatrixBase<DerivedIn>& y0,
		const double tolrel,
		const double tolabs,
		double& x,
		const Eigen::MatrixBase<DerivedOut>& y,
		double& h,
		int& state)
	{
		static_assert(Size > 0, "Size must be a positive integer");
		static_assert((DerivedIn::RowsAtCompileTime == Size) && (DerivedIn::ColsAtCompileTime == 1), "DerivedIn must have compatible size");
		static_assert((DerivedOut::RowsAtCompileTime == Size) && (DerivedOut::ColsAtCompileTime == 1), "DerivedOut must have compatible size");
		typedef typename Eigen::Matrix<double, Size, 1, Options, Size, 1> Vector;

		constexpr double c2 = 0.2;
		constexpr double c3 = 0.3;
		constexpr double c4 = 0.8;
		constexpr double c5 = 8.0 / 9.0;
		constexpr double a21 = 1.0;
		constexpr double a31 = 1.0 / 4.0;
		constexpr double a32 = 3.0 / 4.0;
		constexpr double a41 = 11.0 / 9.0;
		constexpr double a42 = -14.0 / 3.0;
		constexpr double a43 = 40.0 / 9.0;
		constexpr double a51 = 4843.0 / 1458.0;
		constexpr double a52 = -3170.0 / 243.0;
		constexpr double a53 = 8056.0 / 729.0;
		constexpr double a54 = -53.0 / 162.0;
		constexpr double a61 = 9017.0 / 3168.0;
		constexpr double a62 = -355.0 / 33.0;
		constexpr double a63 = 46732.0 / 5247.0;
		constexpr double a64 = 49.0 / 176.0;
		constexpr double a65 = -5103.0 / 18656.0;
		constexpr double a71 = 35.0 / 384.0;
		constexpr double a73 = 500.0 / 1113.0;
		constexpr double a74 = 125.0 / 192.0;
		constexpr double a75 = -2187.0 / 6784.0;
		constexpr double a76 = 11.0 / 84.0;
		constexpr double e1 = 71.0 / 57600.0;
		constexpr double e3 = -71.0 / 16695.0;
		constexpr double e4 = 71.0 / 1920.0;
		constexpr double e5 = -17253.0 / 339200.0;
		constexpr double e6 = 22.0 / 525.0;
		constexpr double e7 = -1.0 / 40.0;

		// Maximal number of allowed steps
		constexpr int maxstep = 100000;
		// Step number after which stiffness detection is activated
		constexpr int nstiff = 1000;
		// Machine epsilon for double precision
		constexpr double eps = 2.3e-16;
		constexpr double safe = 0.9;
		constexpr double fac1 = 0.2;
		constexpr double fac2 = 10.0;
		constexpr double beta = 0.04;
		constexpr double expo1 = 0.2 - beta * 0.75;
		constexpr double facc1 = 1.0 / fac1;
		constexpr double facc2 = 1.0 / fac2;
		const double hmax = x1 - x0;

		// Stiffness detection variables
		int iasti = 0;
		int nonsti = 0;
		double hlamb = 0.0;

		bool reject = false;
		double facold = 1e-4;

		x = x0;
		const_cast<Eigen::MatrixBase<DerivedOut>&>(y) = y0;
		na::internal::dopri5_hinit(f, x0, x1, y0, tolrel, tolabs, h);

		int step = 0;
		int naccept = 0;
		int nreject = 0;
		while (true)
		{
			if (step > maxstep)
			{
				state = -2;
				return;
			}
			if (0.01 * h < std::abs(x) * eps)
			{
				state = -3;
				return;
			}
			bool last = false;
			if (x + 1.01 * h > x1)
			{
				h = x1 - x;
				last = true;
			}
			step += 1;

			Vector k1 = f(x, y);
			Vector k2 = f(x + c2 * h, y + (c2 * h) * (a21 * k1));
			Vector k3 = f(x + c3 * h, y + (c3 * h) * (a31 * k1 + a32 * k2));
			Vector k4 = f(x + c4 * h, y + (c4 * h) * (a41 * k1 + a42 * k2 + a43 * k3));
			Vector k5 = f(x + c5 * h, y + (c5 * h) * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4));
			Vector y6 = y + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5);
			Vector k6 = f(x + h, y6);
			Vector y7 = y + h * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6);
			Vector k7 = f(x + h, y7);

			Vector sk = tolabs + tolrel * y.array().abs().max(y7.array().abs());
			double err = (h * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7)).cwiseQuotient(sk).norm() * std::sqrt(1.0 / Size);
			double hnew = h / std::max<double>(facc2, std::min<double>(facc1, std::pow(err, expo1) / std::pow(facold, beta) / safe));

			if (err < 1.0)
			{
				// Current step is accepted.
				naccept += 1;
				facold = std::max<double>(err, 1e-4);
				// Stiffness detection
				if ((naccept % nstiff == 0) || (iasti > 0))
				{
					double num = (k7 - k6).squaredNorm();
					double den = (y7 - y6).squaredNorm();
					if (den > 0.0)
					{
						hlamb = h * std::sqrt(num / den);
					}
					if (hlamb > 3.25)
					{
						nonsti = 0;
						iasti += 1;
						if (iasti == 15)
						{
							state = -4;
							return;
						}
					}
					else
					{
						nonsti += 1;
						if (nonsti == 6)
						{
							iasti = 0;
						}
					}
				}
				x += h;
				const_cast<Eigen::MatrixBase<DerivedOut>&>(y) = y7;
				if (last)
				{
					h = hnew;
					state = 1;
					return;
				}
				hnew = std::min<double>(hnew, hmax);
				if (reject)
				{
					hnew = std::min<double>(hnew, h);
					reject = false;
				}
			}
			else
			{
				// Current step is rejected.
				hnew = h / std::min<double>(facc1, std::pow(err, expo1) / safe);
				reject = true;
				if (naccept >= 1)
				{
					nreject += 1;
				}
			}
			h = hnew;
		}
	}
}

#endif // !NA_RUNGE_KUTTA_H