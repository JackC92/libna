#include "na/spline/nrb_evaluate.h"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "na/special/gamma.h"
#include "na/spline/bsp_evaluate.h"
#include "na/spline/find_span.h"

namespace na
{
	namespace nurbs
	{
		Eigen::Vector3d evaluate(
			const NURBSEntity& nrb,
			const double u)
		{
			Eigen::Index deg = nrb.degree()[0];
			Eigen::Index span = na::spline::find_span(nrb.knots()[0], deg, u);
			return (nrb.coefs().middleRows(span - deg, deg + 1).transpose() * na::bspline::evaluate_basis(nrb.knots()[0], deg, 0, span, u)).hnormalized();
		}

		Eigen::MatrixX3d evaluate(
			const NURBSEntity& nrb,
			const Eigen::Index m,
			const double u)
		{
			Eigen::Index deg = nrb.degree()[0];
			Eigen::Index span = na::spline::find_span(nrb.knots()[0], deg, u);
			Eigen::Matrix4Xd Cw(4, m + 1);
			for (Eigen::Index k = 0; k <= m; ++k)
			{
				Cw.col(k) = nrb.coefs().middleRows(span - deg, deg + 1).transpose() * na::bspline::evaluate_basis(nrb.knots()[0], deg, k, span, u);
			}
			Eigen::Matrix3Xd Ck(3, m + 1);
			for (Eigen::Index k = 0; k <= m; ++k)
			{
				Ck.col(k) = Cw.block<3, 1>(0, k);
				for (Eigen::Index i = 1; i <= k; ++i)
				{
					Ck.col(k) -= na::gamma::binomial_coefficient(k, i) * Cw.coeffRef(3, i) * Ck.col(k - i);
				}
				Ck.col(k) /= Cw.coeffRef(3, 0);
			}
			return Ck.transpose();
		}

		Eigen::Vector3d evaluate(
			const NURBSEntity& nrb,
			const double u,
			const double v)
		{
			Eigen::Index udeg = nrb.degree()[0];
			Eigen::Index vdeg = nrb.degree()[1];
			Eigen::Index uspan = na::spline::find_span(nrb.knots()[0], udeg, u);
			Eigen::Index vspan = na::spline::find_span(nrb.knots()[1], vdeg, v);
			Eigen::VectorXd Nu = na::bspline::evaluate_basis(nrb.knots()[0], udeg, 0, uspan, u);
			Eigen::VectorXd Nv = na::bspline::evaluate_basis(nrb.knots()[1], vdeg, 0, vspan, v);
			Eigen::Vector4d Sw = Eigen::Vector4d::Zero();
			for (Eigen::Index l = 0; l <= vdeg; ++l)
			{
				Sw += nrb.coefs().middleRows(nrb.number()[0] * (vspan - vdeg + l) + (uspan - udeg), udeg + 1).transpose() * Nu * Nv.coeffRef(l);
			}
			return Sw.hnormalized();
		}

		Eigen::MatrixX3d evaluate(
			const NURBSEntity& nrb,
			const Eigen::Index m,
			const double u,
			const double v)
		{
			Eigen::Index udeg = nrb.degree()[0];
			Eigen::Index vdeg = nrb.degree()[1];
			Eigen::Index uspan = na::spline::find_span(nrb.knots()[0], udeg, u);
			Eigen::Index vspan = na::spline::find_span(nrb.knots()[1], vdeg, v);
			Eigen::MatrixXd Nu(udeg + 1, m + 1);
			Eigen::MatrixXd Nv(vdeg + 1, m + 1);
			Eigen::Matrix4Xd Sw = Eigen::Matrix4Xd::Zero(4, (m + 1) * (m + 2) / 2);
			for (Eigen::Index k = 0; k <= m; ++k)
			{
				Nu.col(k) = na::bspline::evaluate_basis(nrb.knots()[0], udeg, k, uspan, u);
				Nv.col(k) = na::bspline::evaluate_basis(nrb.knots()[1], vdeg, k, vspan, v);
			}
			Eigen::Index index = 0;
			for (Eigen::Index k = 0; k <= m; ++k)
			{
				for (Eigen::Index l = 0; l <= m - k; ++l)
				{
					for (Eigen::Index i = 0; i <= vdeg; ++i)
					{
						Sw.col(index) += nrb.coefs().middleRows(nrb.number()[0] * (vspan - vdeg + l) + (uspan - udeg), udeg + 1).transpose() * Nu.col(k) * Nv.col(l).coeffRef(i);
					}
					index += 1;
				}
			}
			auto offset = [m](const Eigen::Index k, const Eigen::Index l) -> Eigen::Index
			{
				// Sum up (m + 1), m, ..., m + 1 - (k - 1) and l
				return ((m + 1) + (m + 1 - (k - 1))) * k / 2 + l;
			};
			index = 0;
			Eigen::Matrix3Xd Skl(3, (m + 1) * (m + 2) / 2);
			for (Eigen::Index k = 0; k <= m; ++k)
			{
				for (Eigen::Index l = 0; l <= m - k; ++l)
				{
					Skl.col(index) = Sw.block<3, 1>(0, index);
					for (Eigen::Index j = 1; j <= l; ++j)
					{
						Skl.col(index) -= na::gamma::binomial_coefficient(l, j) * Sw.coeffRef(3, j) * Skl.col(index - j);
					}
					for (Eigen::Index i = 1; i <= k; ++i)
					{
						Skl.col(index) -= na::gamma::binomial_coefficient(k, i) * Sw.coeffRef(3, offset(i, 0)) * Skl.col(offset(k - i, l));
						Eigen::Vector3d v2 = Eigen::Vector3d::Zero();
						for (Eigen::Index j = 1; j <= l; ++j)
						{
							v2 += na::gamma::binomial_coefficient(l, j) * Sw.coeffRef(3, offset(i, j)) * Skl.col(offset(k - i, l - j));
						}
						Skl.col(index) -= na::gamma::binomial_coefficient(k, i) * v2;
					}
					Skl.col(index) /= Sw.coeffRef(3, 0);
					index += 1;
				}
			}
			return Skl.transpose();
		}
	}
}
