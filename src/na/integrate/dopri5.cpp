#include "na/integrate/dopri5.h"
#include <cassert>
#include <cmath>
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
		double safe,
		double facl,
		double facr,
		double beta,
		double hmax,
		double h,
		int nmax,
		int nstiff)
	{
		nfcn = 0;
		nstep = 0;
		naccpt = 0;
		nrejct = 0;

		if (xend == x)
		{
			status = 1;
			return;
		}

		// The maximal number of steps
		if (nmax == 0)
		{
			nmax = 100000;
		}
		else if (nmax < 0)
		{
			status = -1;
			return;
		}
		// The step number after which test for stiffness is activated
		if (nstiff == 0)
		{
			nstiff = 1000;
		}
		else if (nstiff < 0)
		{
			// The stiffness test will never be activated
			nstiff = nmax + 10;
		}
		constexpr double uround = std::numeric_limits<double>::epsilon();
		if (safe == 0.0)
		{
			safe = 0.9;
		}
		else if ((safe >= 1.0) || (safe <= 1e-4))
		{
			status = -1;
			return;
		}
		if (facl == 0.0)
		{
			facl = 0.2;
		}
		if (facr == 0.0)
		{
			facr = 10.0;
		}
		if (beta == 0.0)
		{
			beta = 0.04;
		}
		else if (beta < 0.0)
		{
			beta = 0.0;
		}
		else if (beta > 0.2)
		{
			status = -1;
			return;
		}
		if (hmax == 0.0)
		{
			hmax = xend - x;
		}

		constexpr double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0 / 9.0;
		constexpr double a21 = 1.0;
		constexpr double a31 = 0.25, a32 = 0.75;
		constexpr double a41 = 11.0 / 9.0, a42 = -14.0 / 3.0, a43 = 40.0 / 9.0;
		constexpr double a51 = 4843.0 / 1458.0, a52 = -3170.0 / 243.0, a53 = 8056.0 / 729.0, a54 = -53.0 / 162.0;
		constexpr double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;
		constexpr double a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0, a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0;
		constexpr double e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0, e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;

		Eigen::ArrayXd yy1(n), k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), ysti(n);
		Eigen::Map<Eigen::ArrayXd> ymap(y, n);
		Eigen::Map<const Eigen::ArrayXd> rtolmap(rtol, itol ? n : 1), atolmap(atol, itol ? n : 1);

		double posneg = (xend - x > 0.0) ? 1.0 : -1.0;
		double facold = 1e-4;
		double expo1 = 0.2 - beta * 0.75;
		double facc1 = 1.0 / facl;
		double facc2 = 1.0 / facr;

		double atoli = atol[0];
		double rtoli = rtol[0];

		bool last = false;
		double hlamb = 0.0;
		int iasti = 0;

		f(x, y, k1.data());

		hmax = std::abs(hmax);
		if (h == 0.0)
		{
			double dnf;
			double dny;
			if (itol)
			{
				dnf = k1.cwiseQuotient(atolmap + rtolmap.cwiseProduct(ymap.abs())).square().sum();
				dny = ymap.cwiseQuotient(atolmap + rtolmap.cwiseProduct(ymap.abs())).square().sum();
			}
			else
			{
				dnf = k1.cwiseQuotient(atoli + rtoli * ymap.abs()).square().sum();
				dny = ymap.cwiseQuotient(atoli + rtoli * ymap.abs()).square().sum();
			}
			if ((dnf < 1e-10) || (dny <= 1e-10))
			{
				h = 1e-6;
			}
			else
			{
				h = 0.01 * std::sqrt(dny / dnf);
			}
			h = posneg * std::min<double>(h, hmax);

			// Perform an explicit Euler step
			k3 = ymap + h * k1;
			f(x + h, k3.data(), k2.data());

			// Estimate the second derivative of the solution
			double der2;
			if (itol)
			{
				der2 = (k2 - k1).cwiseQuotient(atolmap + rtolmap.cwiseProduct(ymap.abs())).square().sum();
			}
			else
			{
				der2 = (k2 - k1).cwiseQuotient(atoli + rtoli * ymap.abs()).square().sum();
			}
			der2 = std::sqrt(der2) / h;

			// The step size is computed such that
			//     h**5 * max(norm(k1), norm(der2)) = 0.01
			double der12 = std::max<double>(std::abs(der2), std::sqrt(dnf));
			double h1;
			if (der12 <= 1e-15)
			{
				h1 = std::max<double>(1e-6, std::abs(h) * 1e-3);
			}
			else
			{
				h1 = std::pow(0.01 / der12, 0.2);
			}
			h = posneg * std::min<double>({ 100.0 * h, h1, hmax });
		}

		nfcn += 2;
		bool reject = false;
		int nonsti = 0;

		double hold;

		while (true)
		{
			if (nstep > nmax)
			{
				hold = h;
				status = -2;
				return;
			}
			if (0.1 * std::abs(h) <= std::abs(x) * uround)
			{
				hold = h;
				status = -3;
				return;
			}
			if ((x + 1.01 * h - xend) * posneg > 0.0)
			{
				h = xend - x;
				last = true;
			}
			nstep += 1;

			yy1 = ymap + (c2 * h) * (a21 * k1);
			f(x + c2 * h, yy1.data(), k2.data());
			yy1 = ymap + (c3 * h) * (a31 * k1 + a32 * k2);
			f(x + c3 * h, yy1.data(), k3.data());
			yy1 = ymap + (c4 * h) * (a41 * k1 + a42 * k2 + a43 * k3);
			f(x + c4 * h, yy1.data(), k4.data());
			yy1 = ymap + (c5 * h) * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4);
			f(x + c5 * h, yy1.data(), k5.data());
			ysti = ymap + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5);
			f(x + h, ysti.data(), k6.data());
			yy1 = ymap + h * (a71 * k1 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6);
			f(x + h, yy1.data(), k2.data());
			k4 = h * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k2);
			nfcn += 6;

			double err;
			if (itol)
			{
				err = k4.cwiseQuotient(atolmap + rtolmap.cwiseProduct(ymap.abs().cwiseMax(yy1.abs()))).square().sum();
			}
			else
			{
				err = k4.cwiseQuotient(atoli + rtoli * ymap.abs().cwiseMax(yy1.abs())).square().sum();
			}
			err = std::sqrt(err / static_cast<double>(n));

			/// Compute the optimal step size for the next time step
			double fac11 = std::pow(err, expo1);
			// Lund-stabilization
			double fac = fac11 / std::pow(facold, beta);
			// Require facl <= hnew / h <= facr
			fac = std::max<double>(facc2, std::min<double>(facc1, fac / safe));
			double hnew = h / fac;

			if (err <= 1.0)
			{
				// Step accepted
				facold = std::max<double>(err, 1e-4);
				naccpt += 1;

				// Stiffness detection
				if (!(naccpt % nstiff) || (iasti > 0))
				{
					double stnum = (k2 - k6).square().sum();
					double stden = (yy1 - ysti).square().sum();
					if (stden > 0.0)
					{
						hlamb = h * std::sqrt(stnum / stden);
					}
					if (hlamb > 3.25)
					{
						nonsti = 0;
						iasti += 1;
						if (iasti == 15)
						{
							// The problem seems to become stiff at current x
							hold = h;
							status = -4;
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

				k1 = k2;
				ymap = yy1;
				x += h;

				if (last)
				{
					hold = hnew;
					status = 1;
					return;
				}
				if (std::abs(hnew) > hmax)
				{
					hnew = posneg * hmax;
				}
				if (reject)
				{
					hnew = posneg * std::min<double>(std::abs(hnew), std::abs(h));
				}
				reject = false;
			}
			else
			{
				// Step rejected
				reject = true;
				hnew = h / std::min<double>(facc1, fac11 / safe);
				if (naccpt >= 1)
				{
					nrejct += 1;
				}
				last = false;
			}
			h = hnew;
		}
	}

	bool dopri5(
		const std::function<double(double, double)> f,
		const double t0,
		const double t1,
		const double y0,
		const double tolrel,
		const double tolabs,
		double& y)
	{
		int status = 0, nfcn = 0, nstep = 0, naccpt = 0, nrejct = 0;
		double x = t0, xend = t1;
		dopri5(1, [&](const double t, const double* y, double* dydt)
			{
				*dydt = f(t, *y);
			}, x, xend, &y, &tolrel, &tolabs, false, status, nfcn, nstep, naccpt, nrejct);
		return status == 1;
	}

	bool dopri5(
		const std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> f,
		const double t0,
		const double t1,
		const Eigen::VectorXd& y0,
		const double tolrel,
		const double tolabs,
		Eigen::VectorXd& y)
	{
		y = y0;
		int status = 0, nfcn = 0, nstep = 0, naccpt = 0, nrejct = 0;
		double x = t0, xend = t1;
		dopri5(y0.size(), [&](const double t, const double* y, double* dydt)
			{
				Eigen::Map<Eigen::VectorXd>(dydt, y0.size()) = f(t, Eigen::Map<const Eigen::VectorXd>(y, y0.size()));
			}, x, xend, y.data(), &tolrel, &tolabs, false, status, nfcn, nstep, naccpt, nrejct);
		return status == 1;
	}

	bool dopri5(
		const std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> f,
		const double t0,
		const double t1,
		const Eigen::VectorXd& y0,
		const double tolrel,
		const double tolabs,
		Eigen::Ref<Eigen::VectorXd, Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>> y)
	{
		assert((y.size() == y0.size()) && "dopri5: y and y0 must have the same size");

		int status = 0, nfcn = 0, nstep = 0, naccpt = 0, nrejct = 0;
		double x = t0, xend = t1;
		if (y.innerStride() == 1)
		{
			y = y0;
			dopri5(y0.size(), [&](const double t, const double* y, double* dydt)
				{
					Eigen::Map<Eigen::VectorXd>(dydt, y0.size()) = f(t, Eigen::Map<const Eigen::VectorXd>(y, y0.size()));
				}, x, xend, y.data(), &tolrel, &tolabs, false, status, nfcn, nstep, naccpt, nrejct);
		}
		else
		{
			Eigen::VectorXd ytmp = y0;
			dopri5(y0.size(), [&](const double t, const double* y, double* dydt)
				{
					Eigen::Map<Eigen::VectorXd>(dydt, y0.size()) = f(t, Eigen::Map<const Eigen::VectorXd>(y, y0.size()));
				}, x, xend, ytmp.data(), &tolrel, &tolabs, false, status, nfcn, nstep, naccpt, nrejct);
			y = ytmp;
		}
		return status == 1;
	}
}
