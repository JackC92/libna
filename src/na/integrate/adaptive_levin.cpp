#include "na/integrate/adaptive_levin.h"
#include <cmath>
#include <complex>
#include <functional>
#include <queue>
#include <utility>
#include "Eigen/Core"
#include "na/linalg/qr.h"
#include "na/special/chebyshev.h"

namespace na
{
	namespace internal
	{
		template <typename Scalar>
		struct IntegralTriplet
		{
			double a;
			double b;
			Scalar val;
		};

		inline double adaptive_levin_estimate(
			const std::function<double(double)> f,
			const std::function<double(double)> g,
			const double a,
			const double b,
			const Eigen::Index k,
			const Eigen::VectorXd& xscheb,
			const Eigen::MatrixXd& D)
		{
			Eigen::VectorXd x(k), gvec(k), fvec(k), pvec(k);
			Eigen::MatrixXd Q, L;
			x = 0.5 * ((b - a) * xscheb.array() + (b + a));
			for (int i = 0; i < k; ++i)
			{
				fvec(i) = 0.5 * (b - a) * f(x(i));
				gvec(i) = g(x(i));
			}
			na::linalg::qr_factorize(D + (D * gvec).asDiagonal().toDenseMatrix(), Q, L);
			na::linalg::qr_solve(Q, L, fvec, pvec);
			return pvec(k - 1) * std::exp(gvec(k - 1)) - pvec(0) * std::exp(gvec(0));
		}

		inline std::complex<double> cadaptive_levin_estimate(
			const std::function<double(double)> f,
			const std::function<double(double)> g,
			const double a,
			const double b,
			const Eigen::Index k,
			const Eigen::VectorXd& xscheb,
			const Eigen::MatrixXd& D)
		{
			Eigen::VectorXd x(k), gvec(k), fvec(k);
			Eigen::VectorXcd pvec(k);
			Eigen::MatrixXcd Q, L;
			x = 0.5 * ((b - a) * xscheb.array() + (b + a));
			for (int i = 0; i < k; ++i)
			{
				fvec(i) = 0.5 * (b - a) * f(x(i));
				gvec(i) = g(x(i));
			}
			na::linalg::cqr_factorize(D + (std::complex<double>(0.0, 1.0) * (D * gvec)).asDiagonal().toDenseMatrix(), Q, L);
			na::linalg::cqr_solve(Q, L, fvec, pvec);
			return pvec(k - 1) * std::polar(1.0, gvec(k - 1)) - pvec(0) * std::polar(1.0, gvec(0));
		}
	}

	double adaptive_levin(
		const std::function<double(double)> f,
		const std::function<double(double)> g,
		const double a,
		const double b,
		const double tol,
		const Eigen::Index k)
	{
		Eigen::VectorXd xscheb, whtscheb;
		na::chebyshev_pract::quadrature(k, xscheb, whtscheb);
		Eigen::MatrixXd D = na::chebyshev_pract::differentiation_matrix(k, xscheb);

		std::queue<na::internal::IntegralTriplet<double>> q;
		q.push({ a, b, na::internal::adaptive_levin_estimate(f, g, a, b, k, xscheb, D) });
		double val = 0.0;
		while (!q.empty())
		{
			double mid = 0.5 * (q.front().a + q.front().b);
			double val0 = q.front().val;
			double valL = na::internal::adaptive_levin_estimate(f, g, q.front().a, mid, k, xscheb, D);
			double valR = na::internal::adaptive_levin_estimate(f, g, mid, q.front().b, k, xscheb, D);
			if (std::abs(val0 - valL - valR) < tol)
			{
				val += val0;
			}
			else
			{
				q.push({ q.front().a, mid, valL });
				q.push({ mid, q.front().b, valR });
			}
			q.pop();
		}
		return val;
	}

	std::complex<double> cadaptive_levin(
		const std::function<double(double)> f,
		const std::function<double(double)> g,
		const double a,
		const double b,
		const double tol,
		const Eigen::Index k)
	{
		Eigen::VectorXd xscheb, whtscheb;
		na::chebyshev_pract::quadrature(k, xscheb, whtscheb);
		Eigen::MatrixXd D = na::chebyshev_pract::differentiation_matrix(k, xscheb);

		std::queue<na::internal::IntegralTriplet<std::complex<double>>> q;
		q.push({ a, b, na::internal::cadaptive_levin_estimate(f, g, a, b, k, xscheb, D) });
		std::complex<double> val(0.0, 0.0);
		while (!q.empty())
		{
			double mid = 0.5 * (q.front().a + q.front().b);
			std::complex<double> val0 = q.front().val;
			std::complex<double> valL = na::internal::cadaptive_levin_estimate(f, g, q.front().a, mid, k, xscheb, D);
			std::complex<double> valR = na::internal::cadaptive_levin_estimate(f, g, mid, q.front().b, k, xscheb, D);
			if (std::abs(val0 - valL - valR) < tol)
			{
				val += val0;
			}
			else
			{
				q.push({ q.front().a, mid, valL });
				q.push({ mid, q.front().b, valR });
			}
			q.pop();
		}
		return val;
	}
}
