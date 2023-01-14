#include "na/spline/bsp_degelevate.h"
#include <algorithm>
#include "Eigen/Core"
#include "na/core/array.h"
#include "na/special/gamma.h"

namespace na
{
	namespace bspline
	{
		void bsp_degelevate(
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Index deg,
			const Eigen::Index t,
			Eigen::MatrixXd& coefsout,
			Eigen::ArrayXd& knotsout)
		{
			if (t == 0)
			{
				coefsout = coefs;
				knotsout = knots;
				return;
			}

			knotsout = knots;
			na::deduplicate<double>(knotsout);
			coefsout.setZero(coefs.rows() + t * (knotsout.size() - 1), coefs.cols());
			knotsout.setZero(knots.size() + t * knotsout.size());

			// Coefficients for degree elevating the Bezier segments
			Eigen::MatrixXd bezalfs(deg + 1, deg + t + 1);
			// (deg)th-degree Bezier control points of the current segment
			Eigen::MatrixXd bpts(deg + 1, coefs.cols());
			// (deg + t)th-degree Bezier control points of the current segment
			Eigen::MatrixXd ebpts(deg + t + 1, coefs.cols());
			// leftmost control points of the next Bezier segment
			Eigen::MatrixXd nextbpts(deg - 1, coefs.cols());
			// knot insertion alphas
			Eigen::ArrayXd alphas(deg - 1);

			Eigen::Index m = knots.size() - 1;
			Eigen::Index ph = deg + t;
			// Compute Bezier degree elevation coefficients
			// Only the filled in entries are used, so we don't need to initialize the matrix to zero
			Eigen::Index ph2 = ph / 2;
			bezalfs.coeffRef(0, 0) = 1.0;
			bezalfs.coeffRef(deg, ph) = 1.0;
			for (Eigen::Index i = 1; i <= ph2; ++i)
			{
				double inv = 1.0 / na::gamma::binomial_coefficient(ph, i);
				Eigen::Index mpi = std::min<Eigen::Index>(deg, i);
				for (Eigen::Index j = std::max<Eigen::Index>(0, i - t); j <= mpi; ++j)
				{
					// The last term can be particularly large, so the function binomial_coefficient() may need to return a double if we want to allow large t.
					bezalfs.coeffRef(j, i) = inv * na::gamma::binomial_coefficient(deg, j) * na::gamma::binomial_coefficient(t, i - j);
				}
			}
			for (Eigen::Index i = ph2 + 1; i <= ph - 1; ++i)
			{
				Eigen::Index mpi = std::min<Eigen::Index>(deg, i);
				for (Eigen::Index j = std::max<Eigen::Index>(0, i - t); j <= mpi; ++j)
				{
					bezalfs.coeffRef(j, i) = bezalfs.coeffRef(deg - j, ph - i);
				}
			}

			double ua = knots.coeffRef(0);
			double ub = knots.coeffRef(0);
			knotsout.head(ph + 1).setConstant(ua);
			coefsout.row(0) = coefs.row(0);
			Eigen::Index kind = ph + 1;
			Eigen::Index cind = 1;

			// Initialize first Bezier segment
			bpts = coefs.topRows(deg + 1);
			Eigen::Index lbz = 1;
			Eigen::Index rbz = 1;
			Eigen::Index r = -1;
			Eigen::Index a = deg;
			Eigen::Index b = deg + 1;
			while (b < m)
			{
				Eigen::Index i = b;
				Eigen::Index j = b;
				while ((b < m) && (knots.coeffRef(b) == knots.coeffRef(b + 1)))
				{
					b += 1;
				}
				ub = knots.coeffRef(b);
				// a and b are now the index of the last occurence of ua and ub
				Eigen::Index mul = b - i + 1;
				Eigen::Index oldr = r;
				r = deg - mul;
				lbz = (oldr > 0) ? (oldr + 2) / 2 : 1;
				rbz = (r > 0) ? ph - (r + 1) / 2 : ph;
				if (r > 0)
				{
					double numer = ub - ua;
					for (Eigen::Index k = deg; k > mul; --k)
					{
						alphas.coeffRef(k - mul - 1) = numer / (knots.coeffRef(a + k) - ua);
					}
					for (Eigen::Index j = 1; j <= r; ++j)
					{
						Eigen::Index save = r - j;
						Eigen::Index s = mul + j;
						for (Eigen::Index k = deg; k >= s; --k)
						{
							bpts.row(k) = alphas.coeffRef(k - s) * bpts.row(k) + (1.0 - alphas.coeffRef(k - s)) * bpts.row(k - 1);
						}
						nextbpts.row(save) = bpts.row(deg);
					}
				}
				for (Eigen::Index i = lbz; i <= ph; ++i)
				{
					ebpts.row(i).setZero();
					Eigen::Index mpi = std::min<Eigen::Index>(deg, i);
					for (Eigen::Index j = std::max<Eigen::Index>(0, i - t); j <= mpi; ++j)
					{
						ebpts.row(i) += bezalfs.coeffRef(j, i) * bpts.row(j);
					}
				}
				if (oldr > 1)
				{
					Eigen::Index first = kind - 2;
					Eigen::Index last = kind;
					double den = ub - ua;
					double bet = (ub - knotsout.coeffRef(kind - 1)) / den;
					for (Eigen::Index tr = 1; tr < oldr; ++tr)
					{
						i = first;
						j = last;
						Eigen::Index kj = j - kind + 1;
						while (j - i > tr)
						{
							if (i < cind)
							{
								double alf = (ub - knotsout.coeffRef(i)) / (ua - knotsout.coeffRef(i));
								coefsout.row(i) = alf * coefsout.row(i) + (1.0 - alf) * coefsout.row(i - 1);
							}
							if (j >= lbz)
							{
								if (j - tr <= kind - ph + oldr)
								{
									double gam = (ub - knotsout.coeffRef(j - tr)) / den;
									ebpts.row(kj) = gam * ebpts.row(kj) + (1.0 - gam) * ebpts.row(kj + 1);
								}
								else
								{
									ebpts.row(kj) = bet * ebpts.row(kj) + (1.0 - bet) * ebpts.row(kj + 1);
								}
							}
							i += 1;
							j -= 1;
							kj -= 1;
						}
						first -= 1;
						last += 1;
					}
				}
				if (a != deg)
				{
					knotsout.segment(kind, ph - oldr).setConstant(ua);
					kind += ph - oldr;
				}
				coefsout.middleRows(cind, rbz - lbz + 1) = ebpts.middleRows(lbz, rbz - lbz + 1);
				cind += rbz - lbz + 1;
				if (b < m)
				{
					lbz = 1;
					bpts.topRows(r) = nextbpts.topRows(r);
					bpts.bottomRows(deg - r + 1) = coefs.middleRows(b - deg + r, deg - r + 1);
					a = b;
					b = b + 1;
					ua = ub;
				}
				else
				{
					knotsout.segment(kind, ph + 1) = ub;
				}
			}
		}
	}
}
