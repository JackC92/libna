#include "na/linalg/hessenberg.h"
#include <complex>
#include "Eigen/Core"

namespace na
{
	namespace internal
	{
		namespace linalg
		{
			namespace hessenberg
			{
				inline void herm_p_rank1_two_elems_rotate(
					const Eigen::Ref<Eigen::Matrix2cd>& u,
					Eigen::dcomplex& x,
					Eigen::dcomplex& y)
				{
					Eigen::dcomplex d = u(0, 0) * x + u(0, 1) * y;
					y = u(1, 0) * x + u(1, 1) * y;
					x = d;
				}

				inline void herm_p_rank1_two_elems_rotate01(
					const Eigen::Ref<Eigen::Matrix2cd>& u,
					Eigen::dcomplex& x,
					Eigen::dcomplex& y)
				{
					y = u(1, 0) * x + u(1, 1) * y;
				}

				inline void herm_p_rank1_two_elems_rotate10(
					const Eigen::Ref<Eigen::Matrix2cd>& u,
					Eigen::dcomplex& x,
					Eigen::dcomplex& y)
				{
					x = u(0, 0) * x + u(0, 1) * y;
				}

				inline void herm_p_rank1_rotfnd4(
					const Eigen::Vector2cd& a,
					Eigen::Ref<Eigen::Matrix2cd> u)
				{
					const double d = a.norm();
					if (d == 0.0)
					{
						u.setIdentity();
					}
					else
					{
						const double dinv = 1.0 / d;
						u(0, 0) = a(1) * dinv;
						u(1, 1) = std::conj(u(0, 0));
						u(0, 1) = -a(0) * dinv;
						u(1, 0) = -std::conj(u(0, 1));
					}
				}

				inline void herm_p_rank1_iter(
					const int n,
					Eigen::Ref<Eigen::VectorXcd> p,
					Eigen::Ref<Eigen::VectorXcd> q,
					Eigen::Ref<Eigen::VectorXcd> diag,
					Eigen::Ref<Eigen::VectorXcd> supdiag)
				{
					Eigen::Vector2cd xy;
					Eigen::VectorXcd q2, subdiag;
					Eigen::MatrixXcd uus(2 * n - 2, 2);
					q2 = q;
					subdiag = supdiag.conjugate();
					for (int k = n - 1; k > 0; --k)
					{
						xy(0) = supdiag(k - 1) + p(k - 1) * std::conj(q(k));
						xy(1) = diag(k) + p(k) * std::conj(q(k));
						herm_p_rank1_rotfnd4(xy, uus.block<2, 2>(2 * k - 2, 0));
						if (k > 1)
						{
							Eigen::dcomplex subsub = -q2(k) * std::conj(p(k - 2));
							herm_p_rank1_two_elems_rotate10(uus.block<2, 2>(2 * k - 2, 0), subdiag(k - 2), subsub);
						}
						herm_p_rank1_two_elems_rotate(uus.block<2, 2>(2 * k - 2, 0), diag(k - 1), subdiag(k - 1));
						double dd1 = (supdiag(k - 1) * std::conj(supdiag(k - 1)) + diag(k) * std::conj(diag(k))).real();
						double dd2 = (p(k - 1) * std::conj(q(k)) * std::conj(p(k - 1)) * q(k) + p(k) * std::conj(q(k)) * std::conj(p(k)) * q(k)).real();
						if (dd2 < dd1)
						{
							herm_p_rank1_two_elems_rotate01(uus.block<2, 2>(2 * k - 2, 0), supdiag(k - 1), diag(k));
							herm_p_rank1_two_elems_rotate(uus.block<2, 2>(2 * k - 2, 0), p(k - 1), p(k));
						}
						else
						{
							herm_p_rank1_two_elems_rotate(uus.block<2, 2>(2 * k - 2, 0), supdiag(k - 1), diag(k));
							herm_p_rank1_two_elems_rotate01(uus.block<2, 2>(2 * k - 2, 0), p(k - 1), p(k));
							p(k - 1) = -supdiag(k - 1) / std::conj(q(k));
						}
						herm_p_rank1_two_elems_rotate10(uus.block<2, 2>(2 * k - 2, 0), q2(k - 1), q2(k));
					}
					for (int k = n - 1; k > 0; --k)
					{
						supdiag(k - 1) = -p(k - 1) * std::conj(q(k));
						herm_p_rank1_two_elems_rotate(uus.block<2, 2>(2 * k - 2, 0), q(k - 1), q(k));
					}
					uus = uus.conjugate();
					for (int k = n - 1; k > 0; --k)
					{
						herm_p_rank1_two_elems_rotate(uus.block<2, 2>(2 * k - 2, 0), diag(k - 1), supdiag(k - 1));
						herm_p_rank1_two_elems_rotate01(uus.block<2, 2>(2 * k - 2, 0), subdiag(k - 1), diag(k));
					}
				}

				inline void herm_p_rank1_one_dimen(
					Eigen::Ref<Eigen::VectorXcd> diag,
					Eigen::Ref<Eigen::VectorXcd> supdiag,
					Eigen::Ref<Eigen::VectorXcd> p,
					Eigen::Ref<Eigen::VectorXcd> q,
					const int n,
					const double eps)
				{
					Eigen::dcomplex d1, d2, s1, s2, aa, bb, cc, clam, clam1, clam2, clamsum;
					int ifout = 0;
					clamsum = 0;
					for (int iter = 0; iter < 1000; ++iter)
					{
						d1 = diag(0) + p(0) * std::conj(q(0));
						d2 = diag(1) + p(1) * std::conj(q(1));
						s1 = supdiag(0) + p(0) * std::conj(q(1));
						s2 = std::conj(supdiag(0)) + p(1) * std::conj(q(0));
						aa = 1.0;
						bb = -(d1 + d2);
						cc = d1 * d2 - s1 * s2;
						clam1 = (-bb + std::sqrt(bb * bb - 4.0 * aa * cc)) / (2.0 * aa);
						clam2 = (-bb - std::sqrt(bb * bb - 4.0 * aa * cc)) / (2.0 * aa);
						clam = clam1;
						if (std::abs(clam2 - d1) < std::abs(clam1 - d1))
						{
							clam = clam2;
						}
						clamsum += clam;
						for (int i = 0; i < n; ++i)
						{
							diag(i) -= clam;
						}
						herm_p_rank1_iter(n, p, q, diag, supdiag);
						if (std::abs(supdiag(0) + p(0) * std::conj(q(1))) < eps)
						{
							ifout += 1;
						}
						if (ifout == 2)
						{
							break;
						}
					}
					diag.array() += clamsum;
				}
			}
		}
	}

	namespace linalg
	{
		namespace hessenberg
		{
			void herm_p_rank1(
				Eigen::VectorXcd& diag,
				Eigen::VectorXcd& supdiag,
				Eigen::VectorXcd& p,
				Eigen::VectorXcd& q,
				const int n,
				const double eps)
			{
				na::internal::linalg::hessenberg::herm_p_rank1_iter(n, p, q, diag, supdiag);
				na::internal::linalg::hessenberg::herm_p_rank1_iter(n, p, q, diag, supdiag);
				na::internal::linalg::hessenberg::herm_p_rank1_iter(n, p, q, diag, supdiag);
				for (int i = 0; i < n - 1; ++i)
				{
					const int nn = n - i;
					na::internal::linalg::hessenberg::herm_p_rank1_one_dimen(diag.segment(i, nn), supdiag.segment(i, nn - 1), p.segment(i, nn), q.segment(i, nn), nn, eps);
					diag(i) += p(i) * std::conj(q(i));
				}
				diag(n - 1) += p(n - 1) * std::conj(q(n - 1));
			}
		}
	}
}
