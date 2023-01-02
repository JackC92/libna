#include "na/spline/bsp_evaluate.h"
#include <cassert>
#include <cmath>
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "na/integrate/adaptive_simpson.h"
#include "na/spline/find_span.h"

namespace na
{
	namespace bspline
	{
		Eigen::VectorXd evaluate(
			Eigen::Ref<const Eigen::MatrixXd> coefs,
			Eigen::Ref<const Eigen::VectorXd> knots,
			const Eigen::Index deg,
			const Eigen::Index m,
			const double u)
		{
			Eigen::Index mu = na::spline::find_span(knots, deg, u);
			return coefs.middleRows(mu - deg, deg + 1).transpose() * evaluate_basis(knots, deg, m, u);
		}

		Eigen::VectorXd evaluate_basis(
			Eigen::Ref<const Eigen::VectorXd> knots,
			const Eigen::Index deg,
			const Eigen::Index m,
			const double u)
		{
			// n is the number of control points needed to be compatible with T.
			Eigen::Index n = knots.size() - deg - 1;

			assert((0 <= m) && "evaluate_basis: m must be a non-negative integer");
			assert((n >= deg + 1) && "evaluate_basis: T must contain at least 2 * deg + 2 knots");
			assert(((knots(deg) <= u) && (u <= knots(n))) && "evaluate_basis: u must be between [knots(deg), knots(n)]");

			Eigen::Index mu = na::spline::find_span(knots, deg, u);
			Eigen::VectorXd b = Eigen::VectorXd::Zero(deg + 1);
			if (m <= deg)
			{
				b(deg) = 1.0;

				double w1, w2;
				// B_{0} = 1 is a 1x1 matrix and B_{r}(s) = B_{r - 1}(s) * R_{r}(s).
				// Post-multiply by the matrix R_{r}(s).
				for (Eigen::Index r = 1; r <= deg - m; ++r)
				{
					Eigen::Index k = mu - r + 1;
					w2 = (knots(k + r) - u) / (knots(k + r) - knots(k));
					b(deg - r) = w2 * b(deg - r + 1);
					for (Eigen::Index i = deg - r + 1; i <= deg - 1; ++i)
					{
						k += 1;
						w1 = w2;
						w2 = (knots(k + r) - u) / (knots(k + r) - knots(k));
						b(i) = (1.0 - w1) * b(i) + w2 * b(i + 1);
					}
					b(deg) *= 1.0 - w2;
				}

				// Post-multiply by the matrix r * DR_r, which is independent of s because the entries of R_r(s) are all linear in s.
				for (Eigen::Index r = deg - m + 1; r <= deg; ++r)
				{
					Eigen::Index k = mu - r + 1;
					b(deg - r) = -b(deg - r + 1) / (knots(k + r) - knots(k));
					for (Eigen::Index i = deg - r + 1; i <= deg - 1; ++i)
					{
						b(i) = b(i) / (knots(k + r) - knots(k)) - b(i + 1) / (knots(k + r + 1) - knots(k + 1));
						k += 1;
					}
					b(deg) /= knots(k + r) - knots(k);
					b *= r;
				}
			}
			return b;
		}

		Eigen::MatrixXd evaluate_basis(
			Eigen::Ref<const Eigen::VectorXd> knots,
			const Eigen::Index deg,
			const Eigen::Index m,
			Eigen::Ref<const Eigen::VectorXd> u)
		{
			// n is the number of control points needed to be compatible with T.
			const Eigen::Index n = knots.size() - deg - 1;
			// b.row(i) contains the value of all n basis function evaluated at u(i).
			const Eigen::Index num_pts = u.size();
			Eigen::MatrixXd b = Eigen::MatrixXd::Zero(num_pts, n);
			for (Eigen::Index i = 0; i < num_pts; ++i)
			{
				Eigen::Index mu = na::spline::find_span(knots, deg, u(i));
				b.block(i, mu - deg, 1, deg + 1) = evaluate_basis(knots, deg, m, u(i)).transpose();
			}
			return b;
		}

		double length(
			Eigen::Ref<const Eigen::MatrixXd> coefs, 
			Eigen::Ref<const Eigen::VectorXd> knots,
			const Eigen::Index deg,
			const double u0,
			const double u1,
			const double tol)
		{
			return na::adaptive_simpson<double>([&](const double u) -> double
				{
					return evaluate(coefs, knots, deg, 1, u).norm();
				}, u0, u1, tol, 10);
		}

		Eigen::Vector3d curvature_binormal(
			Eigen::Ref<const Eigen::MatrixXd> coefs,
			Eigen::Ref<const Eigen::VectorXd> knots,
			const Eigen::Index deg,
			const double u)
		{
			assert((coefs.cols() == 3) && "curvature_binormal: coefs must be a matrix of size (n, 3)");
			
			Eigen::Vector3d v1 = evaluate(coefs, knots, deg, 1, u);
			Eigen::Vector3d v2 = evaluate(coefs, knots, deg, 2, u);
			return v1.cross(v2) * std::pow(v1.squaredNorm(), -1.5);
		}

		Eigen::Vector3d scaled_curvature_binormal(
			Eigen::Ref<const Eigen::MatrixXd> coefs,
			Eigen::Ref<const Eigen::VectorXd> knots,
			const Eigen::Index deg,
			const double u)
		{
			assert((coefs.cols() == 3) && "scaled_curvature_binormal: coefs must be a matrix of size (n, 3)");

			Eigen::Vector3d v1 = evaluate(coefs, knots, deg, 1, u).normalized();
			Eigen::Vector3d v2 = evaluate(coefs, knots, deg, 2, u);
			return v1.cross(v2);
		}
	}
}
