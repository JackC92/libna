#include "na/spline/nrb_tessellate.h"
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"
#include "na/spline/nrb_evaluate.h"

namespace na
{
	namespace nurbs
	{
		void nrb_tessellate(
			const NURBSEntity& nrb,
			Eigen::MatrixX3d& V,
			const Eigen::Index npts)
		{
			V.resize(npts, 3);
			for (Eigen::Index i = 0; i < npts; ++i)
			{
				V.row(i) = evaluate(nrb, static_cast<double>(i) / (npts - 1.0));
			}
		}

		void nrb_tessellate(
			const NURBSEntity& nrb,
			Eigen::MatrixX3d& V, 
			Eigen::MatrixX3i& F, 
			const Eigen::Index upts,
			const Eigen::Index vpts)
		{
			V.resize(upts * vpts, 3);
			F.resize((upts - 1) * (vpts - 1) * 2, 3);
			double u0 = nrb.knots()[0].coeffRef(nrb.degree()[0]);
			double u1 = nrb.knots()[0].coeffRef(nrb.number()[0]);
			double v0 = nrb.knots()[1].coeffRef(nrb.degree()[1]);
			double v1 = nrb.knots()[1].coeffRef(nrb.number()[1]);
			for (Eigen::Index j = 0; j < vpts; ++j)
			{
				double v = v0 + static_cast<double>(j) / (vpts - 1.0) * (v1 - v0);
				for (Eigen::Index i = 0; i < upts; ++i)
				{
					double u = u0 + static_cast<double>(i) / (upts - 1.0) * (u1 - u0);
					V.row(j * upts + i) = evaluate(nrb, u, v);
				}
			}
			for (Eigen::Index j = 0; j < vpts - 1; ++j)
			{
				F.block<Eigen::Dynamic, 1>(j * (upts - 1) * 2, 0, upts - 1, 1).setLinSpaced(j * upts + 0, (j + 1) * upts - 2);
				F.block<Eigen::Dynamic, 1>(j * (upts - 1) * 2, 1, upts - 1, 1).setLinSpaced(j * upts + 1, (j + 1) * upts - 1);
				F.block<Eigen::Dynamic, 1>(j * (upts - 1) * 2, 2, upts - 1, 1).setLinSpaced((j + 1) * upts, (j + 2) * upts - 2);
				F.block<Eigen::Dynamic, 1>(j * (upts - 1) * 2 + (upts - 1), 0, upts - 1, 1).setLinSpaced(j * upts + 1, (j + 1) * upts - 1);
				F.block<Eigen::Dynamic, 1>(j * (upts - 1) * 2 + (upts - 1), 1, upts - 1, 1).setLinSpaced((j + 1) * upts + 1, (j + 2) * upts - 1);
				F.block<Eigen::Dynamic, 1>(j * (upts - 1) * 2 + (upts - 1), 2, upts - 1, 1).setLinSpaced((j + 1) * upts + 0, (j + 2) * upts - 2);
			}
		}
	}
}
