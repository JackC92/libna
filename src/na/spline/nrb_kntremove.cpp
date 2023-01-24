#include "na/spline/nrb_kntremove.h"
#include <algorithm>
#include "Eigen/Core"
#include "na/spline/nrb_entity.h"

namespace na
{
	namespace nurbs
	{
		Eigen::Index kntremove_curve(
			NURBSEntity& nrb,
			const Eigen::Index r,
			const Eigen::Index num,
			const double tol)
		{
			// r should be the index of the last occurence of the knot
			double u = nrb.knots()[0].coeffRef(r);
			Eigen::Index m = nrb.knots()[0].size() - 1;
			Eigen::Index n = nrb.number()[0] - 1;
			Eigen::Index s = std::count(nrb.knots()[0].data(), nrb.knots()[0].data() + nrb.knots()[0].size(), nrb.knots()[0].coeffRef(r));
			Eigen::Index deg = nrb.degree()[0];
			Eigen::Index ord = deg + 1;
			Eigen::Index fout = (2 * r - s - deg) / 2;
			Eigen::Index last = r - s;
			Eigen::Index first = r - deg;
			Eigen::MatrixX4d temp(2 * deg + 1, 4);
			Eigen::Index t = 0;
			for (t = 0; t < num; ++t)
			{
				Eigen::Index off = first - 1;
				temp.row(0) = nrb.coefs().row(off);
				temp.row(last + 1 - off) = nrb.coefs().row(last + 1);
				Eigen::Index i = first;
				Eigen::Index j = last;
				Eigen::Index ii = 1;
				Eigen::Index jj = last - off;
				bool remflag = false;
				while (j - i > t)
				{
					double alfi = (u - nrb.knots()[0].coeffRef(i)) / (nrb.knots()[0].coeffRef(i + ord + t) - nrb.knots()[0].coeffRef(i));
					double alfj = (u - nrb.knots()[0].coeffRef(j - t)) / (nrb.knots()[0].coeffRef(j + ord) - nrb.knots()[0].coeffRef(j - t));
					temp.row(ii) = (nrb.coefs().row(i) - (1.0 - alfi) * temp.row(ii - 1)) / alfi;
					temp.row(jj) = (nrb.coefs().row(j) - alfj * temp.row(jj + 1)) / (1.0 - alfj);
					i += 1;
					j -= 1;
					ii += 1;
					jj -= 1;
				}
				if (j - i < t)
				{
					if ((temp.row(ii - 1) - temp.row(jj + 1)).norm() <= tol)
					{
						remflag = true;
					}
				}
				else
				{
					double alfi = (u - nrb.knots()[0].coeffRef(i)) / (nrb.knots()[0].coeffRef(i + ord + t) - nrb.knots()[0].coeffRef(i));
					if ((nrb.coefs().row(i) - (alfi * temp.row(ii + t + 1) + (1.0 - alfi) * temp.row(ii - 1))).norm() <= tol)
					{
						remflag = true;
					}
				}
				// Ignore the tolerance and remove the knot
				if (tol <= 0.0)
				{
					remflag = true;
				}
				if (!remflag)
				{
					break;
				}
				else
				{
					i = first;
					j = last;
					while (j - i > t)
					{
						nrb.coefs().row(i) = temp.row(i - off);
						nrb.coefs().row(j) = temp.row(j - off);
						i += 1;
						j -= 1;
					}
				}
				first -= 1;
				last += 1;
			}
			if (t == 0)
			{
				return t;
			}
			for (Eigen::Index k = r + 1; k <= m; ++k)
			{
				nrb.knots()[0].coeffRef(k - t) = nrb.knots()[0].coeffRef(k);
			}
			Eigen::Index j = fout;
			Eigen::Index i = j;
			for (Eigen::Index k = 1; k < t; ++k)
			{
				if ((k % 2) == 1)
				{
					i += 1;
				}
				else
				{
					j -= 1;
				}
			}
			for (Eigen::Index k = i + 1; k <= n; ++k)
			{
				nrb.coefs().row(j) = nrb.coefs().row(k);
				j += 1;
			}
			nrb.coefs().conservativeResize(nrb.coefs().rows() - t, nrb.coefs().cols());
			nrb.knots()[0].conservativeResize(m + 1 - t);
			nrb.number()[0] = nrb.coefs().rows();
			return t;
		}
	}
}
