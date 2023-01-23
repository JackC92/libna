#include "na/spline/bsp_evaluate.h"
#include <cassert>
#include <cmath>
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/Sparse"
#include "na/integrate/adaptive_simpson.h"
#include "na/spline/find_span.h"

namespace na
{
	namespace bspline
	{
		Eigen::VectorXd evaluate(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index deg,
			const Eigen::Index m,
			const double u)
		{
			Eigen::Index span = na::spline::find_span(knots, deg, u);
			return coefs.middleRows(span - deg, deg + 1).transpose() * evaluate_basis(knots, deg, m, span, u);
		}

		Eigen::VectorXd evaluate_basis(
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index deg,
			const Eigen::Index m,
			const Eigen::Index span,
			const double u)
		{
			// n is the number of control points needed to be compatible with T.
			Eigen::Index n = knots.size() - deg - 1;

			assert((0 <= m) && "evaluate_basis: m must be a non-negative integer");
			assert((n >= deg + 1) && "evaluate_basis: T must contain at least 2 * deg + 2 knots");
			assert(((knots.coeffRef(deg) <= u) && (u <= knots.coeffRef(n))) && "evaluate_basis: u must be between [knots(deg), knots(n)]");

			Eigen::VectorXd N = Eigen::VectorXd::Zero(deg + 1);
			if (m <= deg)
			{
				N.coeffRef(deg) = 1.0;

				double w1, w2;
				// B_{0} = 1 is a 1x1 matrix and B_{r}(s) = B_{r - 1}(s) * R_{r}(s).
				// Post-multiply by the matrix R_{r}(s).
				for (Eigen::Index r = 1; r <= deg - m; ++r)
				{
					Eigen::Index k = span - r + 1;
					w2 = (knots.coeffRef(span + 1) - u) / (knots.coeffRef(span + 1) - knots.coeffRef(span - r + 1));
					N.coeffRef(deg - r) = w2 * N.coeffRef(deg - r + 1);
					for (Eigen::Index i = deg - r + 1; i <= deg - 1; ++i)
					{
						k += 1;
						w1 = w2;
						w2 = (knots.coeffRef(k + r) - u) / (knots.coeffRef(k + r) - knots.coeffRef(k));
						N.coeffRef(i) = (1.0 - w1) * N.coeffRef(i) + w2 * N.coeffRef(i + 1);
					}
					N.coeffRef(deg) *= 1.0 - w2;
				}

				// Post-multiply by the matrix r * DR_r, which is independent of s because the entries of R_r(s) are all linear in s.
				w1 = 1.0;
				for (Eigen::Index r = deg - m + 1; r <= deg; ++r)
				{
					N.coeffRef(deg - r) = -N.coeffRef(deg - r + 1) / (knots.coeffRef(span + 1) - knots.coeffRef(span - r + 1));
					N.segment(deg - r + 1, r - 1) = N.segment(deg - r + 1, r - 1).array().cwiseQuotient(knots.segment(span + 1, r - 1) - knots.segment(span - r + 1, r - 1)) - N.segment(deg - r + 2, r - 1).array().cwiseQuotient(knots.segment(span + 2, r - 1) - knots.segment(span - r + 2, r - 1));
					N.coeffRef(deg) /= knots.coeffRef(span + r) - knots.coeffRef(span);
					w1 *= r;
				}
				N *= w1;
			}
			return N;
		}

		Eigen::SparseMatrix<double, Eigen::RowMajor> evaluate_basis(
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index deg,
			const Eigen::Index m,
			const Eigen::Ref<const Eigen::ArrayXi>& span,
			const Eigen::Ref<const Eigen::ArrayXd>& u)
		{
			Eigen::SparseMatrix<double, Eigen::RowMajor> N(u.size(), knots.size() - deg - 1);
			N.data().reserve(u.size() * (deg + 1));
			Eigen::Map<Eigen::ArrayXi, Eigen::Aligned16>(N.outerIndexPtr(), u.size() + 1) = Eigen::ArrayXi::LinSpaced(u.size() + 1, 0, u.size() * (deg + 1));
			for (Eigen::Index k = 0; k < u.size(); ++k)
			{
				Eigen::Map<Eigen::ArrayXi, Eigen::Unaligned>(N.innerIndexPtr() + k * (deg + 1), deg + 1) = Eigen::ArrayXi::LinSpaced(deg + 1, span.coeffRef(k) - deg, span.coeffRef(k));
				Eigen::Map<Eigen::ArrayXd, Eigen::Aligned16>(N.valuePtr() + k * (deg + 1), deg + 1) = evaluate_basis(knots, deg, m, span.coeffRef(k), u.coeffRef(k));
			}
			N.data().resize(u.size() * (deg + 1));
			return N;
		}

		double length(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index deg,
			const double u0,
			const double u1,
			const double tol)
		{
			return na::adaptive_simpson<double>([&](const double u) -> double
				{
					return evaluate(coefs, knots, deg, 1, u).norm();
				}, u0, u1, tol, 20);
		}

		Eigen::Vector3d curvature_binormal(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index deg,
			const double u)
		{
			assert((coefs.cols() == 3) && "curvature_binormal: coefs must be a matrix of size (n, 3)");

			Eigen::Vector3d v1 = evaluate(coefs, knots, deg, 1, u);
			Eigen::Vector3d v2 = evaluate(coefs, knots, deg, 2, u);
			return v1.cross(v2) * std::pow(v1.squaredNorm(), -1.5);
		}

		Eigen::Vector3d scaled_curvature_binormal(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index deg,
			const double u)
		{
			assert((coefs.cols() == 3) && "scaled_curvature_binormal: coefs must be a matrix of size (n, 3)");

			Eigen::Vector3d v1 = evaluate(coefs, knots, deg, 1, u);
			Eigen::Vector3d v2 = evaluate(coefs, knots, deg, 2, u);
			return v1.cross(v2) / v1.squaredNorm();
		}
	}
}
