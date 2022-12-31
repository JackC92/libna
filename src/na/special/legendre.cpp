#include "na/special/legendre.h"
#include <cassert>
#include <cmath>
#include "Eigen/Core"
#include "na/special/bessel.h"

namespace na
{
	namespace detail
	{
		namespace legendre
		{
			inline void evaluate_taylor_series(
				const double x,
				const double h,
				const Eigen::Index n,
				const Eigen::Index k,
				double& pol,
				double& der)
			{
				static bool called = false;
				static Eigen::VectorXd aux(120);

				if (!called)
				{
					for (Eigen::Index i = 0; i < 30; ++i)
					{
						aux.coeffRef(4 * i) = 1.0 / ((i + 2.0) * (i + 3.0));
						aux.coeffRef(4 * i + 1) = (i + 2.0) * (i + 2.0);
						aux.coeffRef(4 * i + 2) = (i + 1.0) * (i + 2.0);
						aux.coeffRef(4 * i + 3) = i + 3.0;
					}
					called = true;
				}

				const double dn = n * (n + 1.0);
				const double d7 = h * h / (1.0 - x * x);
				const double hinv = 1.0 / h;
				const double dd = 2.0 * x * hinv;
				double q1 = der * h;
				double q2 = (x * der - 0.5 * dn * pol) * d7;
				pol = pol + q1 + q2;
				der = (q1 + q2 + q2) * hinv;

				double* ptr = aux.data(), * end = aux.data() + 4 * (k - 2);
				while (ptr != end)
				{
					double d = ptr[0] * (ptr[1] * dd * q2 - (dn - ptr[2]) * q1) * d7;
					der += d * hinv * ptr[3];
					pol += d;
					q1 = q2;
					q2 = d;
					ptr = ptr + 4;
				}
			}

			// This function is based on the paper
			//   A. Glaser, X. Liu, and V. Rokhlin, "A Fast Algorithm for the Calculation of the Roots of Special Functions",
			//   SIAM Journal on Scientific Computing, 2007
			void quadrature_glaser(
				const Eigen::Index n,
				Eigen::VectorXd& xslege,
				Eigen::VectorXd& whtslege)
			{
				if (n == 1)
				{
					xslege.coeffRef(0) = 0.0;
					whtslege.coeffRef(0) = 2.0;
					return;
				}

				const Eigen::Index i = n / 2;
				const Eigen::Index ifodd = n % 2;
				const double n0 = M_PI / (n + 0.5);
				const double n1 = 0.125 * (n - 1) * std::pow(n, -3) + 0.1015625 * std::pow(n, -4) - 1.0;
				const double n2 = 0.07291666666666667 * std::pow(n, -4);

				double* xptr = xslege.data(), * whtptr = whtslege.data();
				for (Eigen::Index j = i + 1; j <= n; ++j)
				{
					double temp = std::cos(n0 * (j - 0.25));
					*xptr++ = temp * (n2 / (temp * temp - 1.0) + n1);
				}
				xptr -= i + ifodd;
				double x0 = 0.0, x1 = xptr[0];
				double pol, der;
				na::legendre::evaluate_polynomial(n, 0.0, pol, der);
				if (ifodd == 1)
				{
					xptr[0] = 0.0;
					whtslege[0] = der;
					x0 = x1;
					x1 = xptr[1];
				}
				for (Eigen::Index j = ifodd; j < i + ifodd; ++j)
				{
					int ifstop = 0;
					for (int k = 0; k < 10; ++k)
					{
						double h = x1 - x0;
						if ((k > 0) && (std::abs(h) < 1e-36))
						{
							break;
						}
						evaluate_taylor_series(x0, h, n, ((k > 0) ? ((k == 1) ? 15 : 6) : 30), pol, der);
						x0 = x1;
						x1 -= pol / der;
						if (std::abs(pol) < 1e-8)
						{
							if (ifstop == 2)
							{
								break;
							}
							ifstop += 1;
						}
					}
					xptr[j] = x1;
					whtptr[j] = der;
					x0 = x1;
					x1 = xptr[j + 1];
				}
				// The non-negative (including zero if the polynomial is odd) nodes and weights have been computed.
				// Move them to the end of the vectors and perform a reverse copy due to symmetry of the nodes and weights.
				for (double* curr = xptr + i + ifodd - 1, *target = xptr + n - 1; curr >= xptr; --curr, --target)
				{
					*target = *curr;
				}
				for (double* curr = whtptr + i + ifodd - 1, *target = whtptr + n - 1; curr >= whtptr; --curr, --target)
				{
					*target = *curr;
				}
				for (double* curr = xptr + n - i, *target = xptr + i - 1; curr < xptr + n; ++curr, --target)
				{
					*target = -(*curr);
				}
				for (double* curr = whtptr + n - i, *target = whtptr + i - 1; curr < whtptr + n; ++curr, --target)
				{
					*target = *curr;
				}
				whtslege = 2.0 / ((1.0 - xslege.cwiseAbs2().array()) * whtslege.cwiseAbs2().array());
			}

			// This function is based on the paper
			//   I. Bogaert, "Iteration-Free Computation of Gauss-Legendre Quadrature Nodes and Weights",
			//   SIAM Journal on Scientific Computing, 2014
			void quadrature_bogaert(
				const Eigen::Index n,
				Eigen::VectorXd& xslege,
				Eigen::VectorXd& whtslege)
			{
				assert((n > 100) && "quadrature_bogaert: n must be greater than 100 to yield machine precision in the IEEE 754 double-precision format");

				const Eigen::Index i = n / 2;
				const Eigen::Index ifodd = n % 2;

				const double vn = 1.0 / (n + 0.5);
				const double vn2 = vn * vn;

				double* xptr = xslege.data(), * whtptr = whtslege.data();
				for (Eigen::Index j = 0; j < i; ++j)
				{
					const double j0k = na::bessel::besselj0_zero(j + 1);
					const double alpha = vn * j0k;
					const double alpha2 = alpha * alpha;
					const double cotalpha = std::cos(alpha) / std::sin(alpha);

					const double g1 = (((((-1.29052996274280508473467968379e-12 * alpha2 + 2.40724685864330121825976175184e-10) * alpha2 - 3.13148654635992041468855740012e-8) * alpha2 + 0.275573168962061235623801563453e-5) * alpha2 - 0.148809523713909147898955880165e-3) * alpha2 + 0.416666666665193394525296923981e-2) * alpha2 - 0.416666666666662959639712457549e-1;
					const double g2 = (((((+2.20639421781871003734786884322e-9 * alpha2 - 7.53036771373769326811030753538e-8) * alpha2 + 0.161969259453836261731700382098e-5) * alpha2 - 0.253300326008232025914059965302e-4) * alpha2 + 0.282116886057560434805998583817e-3) * alpha2 - 0.209022248387852902722635654229e-2) * alpha2 + 0.815972221772932265640401128517e-2;
					const double g3 = (((((-2.97058225375526229899781956673e-8 * alpha2 + 5.55845330223796209655886325712e-7) * alpha2 - 0.567797841356833081642185432056e-5) * alpha2 + 0.418498100329504574443885193835e-4) * alpha2 - 0.251395293283965914823026348764e-3) * alpha2 + 0.128654198542845137196151147483e-2) * alpha2 - 0.416012165620204364833694266818e-2;

					const double h1 = ((((((((-2.20902861044616638398573427475e-14 * alpha2 + 2.30365726860377376873232578871e-12) * alpha2 - 1.75257700735423807659851042318e-10) * alpha2 + 1.03756066927916795821098009353e-8) * alpha2 - 4.63968647553221331251529631098e-7) * alpha2 + 0.149644593625028648361395938176e-4) * alpha2 - 0.326278659594412170300449074873e-3) * alpha2 + 0.436507936507598105249726413120e-2) * alpha2 - 0.305555555555553028279487898503e-1) * alpha2 + 0.833333333333333302184063103900e-1;
					const double h2 = (((((((+3.63117412152654783455929483029e-12 * alpha2 + 7.67643545069893130779501844323e-11) * alpha2 - 7.12912857233642220650643150625e-9) * alpha2 + 2.11483880685947151466370130277e-7) * alpha2 - 0.381817918680045468483009307090e-5) * alpha2 + 0.465969530694968391417927388162e-4) * alpha2 - 0.407297185611335764191683161117e-3) * alpha2 + 0.268959435694729660779984493795e-2) * alpha2 - 0.111111111111214923138249347172e-1;
					const double h3 = (((((((+2.01826791256703301806643264922e-9 * alpha2 - 4.38647122520206649251063212545e-8) * alpha2 + 5.08898347288671653137451093208e-7) * alpha2 - 0.397933316519135275712977531366e-5) * alpha2 + 0.200559326396458326778521795392e-4) * alpha2 - 0.422888059282921161626339411388e-4) * alpha2 - 0.105646050254076140548678457002e-3) * alpha2 - 0.947969308958577323145923317955e-4) * alpha2 + 0.656966489926484797412985260842e-2;

					const double d1 = j0k / std::sin(alpha);
					const double d2 = na::bessel::besselj1_squared(j + 1) * d1;
					const double d3 = vn2 * d1;
					const double d4 = d3 * d3;

					xptr[j] = std::cos(vn * (j0k + alpha * d3 * (g1 + d4 * (g2 + d4 * g3))));
					whtptr[j] = (2.0 * vn) / (d2 + d2 * d4 * (h1 + d4 * (h2 + d4 * h3)));
				}
				for (double* curr = xptr, *target = xptr + n - 1; curr < xptr + i; ++curr, --target)
				{
					*target = *curr;
					*curr = -(*curr);
				}
				for (double* curr = whtptr, *target = whtptr + n - 1; curr < whtptr + i; ++curr, --target)
				{
					*target = *curr;
				}
				if (ifodd == 1)
				{
					xptr[i] = 0.0;
					whtptr[i] = 2.0 / std::pow(na::legendre::evaluate_derivative(n, 0.0), 2);
				}
			}
		}
	}

	namespace legendre
	{
		Eigen::VectorXd polynomial_coefficients(const Eigen::Index n)
		{
			assert((n >= 0) && "polynomial_coefficients: n must be a non-negative integer");

			Eigen::VectorXd coefs(n + 1);
			if (n == 0)
			{
				coefs(0) = 1.0;
			}
			else if (n == 1)
			{
				coefs(0) = 0.0;
				coefs(1) = 1.0;
			}
			else
			{
				// Legendre polynomials can be found using Bonnet's Recursion Formula: n * P_n(x) = (2n - 1) * x * P_{n-1}(x) - (n - 1) * P_{n-2}(x) for n >= 2.
				Eigen::VectorXd coefs1(n + 1), coefs2(n + 1);
				coefs1(0) = 1.0;
				coefs2(0) = 0.0;
				coefs2(1) = 1.0;
				for (Eigen::Index i = 2; i <= n; ++i)
				{
					double a = 2.0 - 1.0 / i, b = 1.0 - 1.0 / i;
					for (Eigen::Index j = 1; j <= i - 2; ++j)
					{
						coefs(j) = a * coefs2(j - 1) - b * coefs1(j);
					}
					coefs(0) = -b * coefs1(0);
					coefs(i) = a * coefs2(i - 1);
					coefs(i - 1) = a * coefs2(i - 2);
					if (i < n)
					{
						std::copy(coefs2.begin(), coefs2.begin() + i, coefs1.begin());
						std::copy(coefs.begin(), coefs.begin() + i + 1, coefs2.begin());
					}
				}
			}
			return coefs;
		}

		double evaluate_polynomial(
			const Eigen::Index n,
			const double x)
		{
			assert((n >= 0) && "evaluate_polynomial: n must be a non-negative integer");

			if (n == 0)
			{
				return 1.0;
			}
			double pkm1, pk, pkp1;
			pk = 1.0;
			pkp1 = x;
			for (Eigen::Index k = 1; k < n; ++k)
			{
				pkm1 = pk;
				pk = pkp1;
				pkp1 = ((2.0 * k + 1.0) * x * pk - k * pkm1) / (k + 1.0);
			}
			return pkp1;
		}

		double evaluate_derivative(
			const Eigen::Index n,
			const double x)
		{
			assert((n >= 0) && "evaluate_derivative: n must be a non-negative integer");

			if (n == 0)
			{
				return 0.0;
			}
			// The derivative is computed using the formula P_{n+1}'(x) = (n+1) * P_{n}(x) + x * P_{n}'(x).
			// Note that this formula holds at the endpoints, but requires O(n) operations.
			double val, pkm1, pk, pkp1;
			val = 1.0;
			pk = 1.0;
			pkp1 = x;
			for (Eigen::Index k = 1; k < n; ++k)
			{
				pkm1 = pk;
				pk = pkp1;
				pkp1 = ((2.0 * k + 1.0) * x * pk - k * pkm1) / (k + 1.0);
				val = (k + 1.0) * pk + x * val;
			}
			return val;
		}

		void evaluate_polynomial(
			const Eigen::Index n,
			const double x,
			double& pol,
			double& der)
		{
			assert((n >= 0) && "evaluate_polynomial: n must be a non-negative integer");

			if (n == 0)
			{
				pol = 1.0;
				der = 0.0;
				return;
			}

			der = 1.0;
			double pkm1, pk, pkp1;
			pk = 1.0;
			pkp1 = x;
			for (Eigen::Index k = 1; k < n; ++k)
			{
				pkm1 = pk;
				pk = pkp1;
				pkp1 = ((2.0 * k + 1.0) * x * pk - k * pkm1) / (k + 1.0);
				der = (k + 1.0) * pk + x * der;
			}
			pol = pkp1;
		}

		void quadrature(
			const Eigen::Index n,
			Eigen::VectorXd& xslege,
			Eigen::VectorXd& whtslege)
		{
			assert((n >= 1) && "quadrature: n must be a positive integer");
			xslege.resize(n);
			whtslege.resize(n);
			if (n > 100)
			{
				na::detail::legendre::quadrature_bogaert(n, xslege, whtslege);
			}
			else
			{
				na::detail::legendre::quadrature_glaser(n, xslege, whtslege);
			}
		}

		Eigen::VectorXd expansion_coefficients(
			const Eigen::Index n,
			const Eigen::VectorXd& vals)
		{
			assert((n >= 1) && "expansion_coefficients: n must be a positive integer");

			Eigen::VectorXd xslege, whtslege;
			quadrature(n, xslege, whtslege);
			return coefficient_matrix(n, xslege, whtslege) * vals;
		}

		double evaluate_expansion(
			const Eigen::Index n,
			const Eigen::VectorXd& coefs,
			const double x)
		{
			assert((n >= 1) && "evaluate_expansion: n must be a positive integer");

			double val = coefs.coeffRef(0) * std::sqrt(0.5);
			if (n > 1)
			{
				val += coefs.coeffRef(1) * std::sqrt(1.5) * x;
				double pkm1, pk, pkp1;
				pk = 1.0;
				pkp1 = x;
				for (Eigen::Index k = 1; k < n - 1; ++k)
				{
					pkm1 = pk;
					pk = pkp1;
					pkp1 = ((2.0 * k + 1.0) * x * pk - k * pkm1) / (k + 1.0);
					val += coefs.coeffRef(k + 1) * std::sqrt(k + 1.5) * pkp1;
				}
			}
			return val;
		}

		void evaluate_expansion(
			const Eigen::Index n,
			const Eigen::VectorXd& coefs,
			const double x,
			double& val,
			double& der)
		{
			assert((n >= 1) && "evaluate_expansion: n must be a positive integer");

			val = coefs.coeffRef(0) * std::sqrt(0.5);
			der = 0.0;
			if (n > 1)
			{
				val += coefs.coeffRef(1) * std::sqrt(1.5) * x;
				der += coefs.coeffRef(1) * std::sqrt(1.5);
				double pkm1, pk, pkp1, temp;
				pk = 1.0;
				pkp1 = x;
				temp = 1.0;
				for (Eigen::Index k = 1; k < n - 1; ++k)
				{
					pkm1 = pk;
					pk = pkp1;
					pkp1 = ((2.0 * k + 1.0) * x * pk - k * pkm1) / (k + 1.0);
					temp = (k + 1.0) * pk + x * temp;
					val += coefs.coeffRef(k + 1) * std::sqrt(k + 1.5) * pkp1;
					der += coefs.coeffRef(k + 1) * std::sqrt(k + 1.5) * temp;
				}
			}
		}

		double interpolate(
			const Eigen::Index n,
			const Eigen::VectorXd& xslege,
			const Eigen::VectorXd& whtslege,
			const Eigen::VectorXd& vals,
			const double x)
		{
			assert((n >= 1) && "interpolate: n must be a positive integer");

			for (Eigen::Index i = 0; i < n; ++i)
			{
				if (x == xslege(i))
				{
					return vals(i) / std::sqrt(whtslege(i));
				}
			}
			Eigen::VectorXd whtsinterp(n);
			for (Eigen::Index i = 0; i < n; ++i)
			{
				// This barycentric weight differs from the literatures by sqrt(w_i) because it's factored into vals.
				// But it will require multiyplying by sqrt(w_i) in the denominator.
				whtsinterp(i) = ((i % 2) ? -1.0 : 1.0) * std::sqrt(1 - std::pow(xslege(i), 2)) / (x - xslege(i));
			}
			return whtsinterp.dot(vals) / whtsinterp.dot(whtslege.cwiseSqrt());
		}

		Eigen::MatrixXd interpolation_matrix(
			const Eigen::Index n,
			const Eigen::VectorXd& xslege,
			const Eigen::VectorXd& whtslege,
			const Eigen::Index m,
			const Eigen::VectorXd& xsout,
			const Eigen::VectorXd& whtsout)
		{
			assert((n >= 1) && "interpolation_matrix: n must be a positive integer");
			assert((m >= 1) && "interpolation_matrix: m must be a positive integer");

			Eigen::MatrixXd ainterp = Eigen::MatrixXd::Zero(m, n);

			Eigen::VectorXd lambda(n);
			for (Eigen::Index j = 0; j < n; ++j)
			{
				lambda(j) = ((j % 2) ? -1.0 : 1.0) * std::sqrt(1.0 - std::pow(xslege(j), 2));
			}
			for (Eigen::Index i = 0; i < m; ++i)
			{
				bool found = false;
				double coeff = 0.0;
				for (int j = 0; j < n; ++j)
				{
					if (xsout(i) == xslege(j))
					{
						ainterp(i, j) = std::sqrt(whtsout(i)) / std::sqrt(whtslege(j));
						found = true;
						break;
					}
					else
					{
						coeff += lambda(j) * std::sqrt(whtslege(j)) / (xsout(i) - xslege(j));
					}
				}
				if (found)
				{
					continue;
				}
				coeff = std::sqrt(whtsout(i)) / coeff;
				for (Eigen::Index j = 0; j < n; ++j)
				{
					ainterp(i, j) = coeff * lambda(j) / (xsout(i) - xslege(j));
				}
			}
			return ainterp;
		}

		Eigen::MatrixXd coefficient_matrix(
			const Eigen::Index n,
			const Eigen::VectorXd& xslege,
			const Eigen::VectorXd& whtslege)
		{
			assert((n >= 1) && "coefficient_matrix: n must be a positive integer");

			Eigen::MatrixXd umatr = Eigen::MatrixXd::Zero(n, n);

			for (Eigen::Index j = 0; j < n; ++j)
			{
				double wht = std::sqrt(whtslege(j));
				umatr(0, j) = wht * std::sqrt(0.5);
				if (n > 1)
				{
					umatr(1, j) = wht * std::sqrt(1.5) * xslege(j);
					double x, pkm1, pk, pkp1;
					x = xslege(j);
					pk = 1.0;
					pkp1 = x;
					for (Eigen::Index k = 1; k < n - 1; ++k)
					{
						pkm1 = pk;
						pk = pkp1;
						pkp1 = ((2.0 * k + 1.0) * x * pk - k * pkm1) / (k + 1.0);
						umatr(k + 1, j) = wht * std::sqrt(k + 1.5) * pkp1;
					}
				}
			}
			return umatr;
		}

		Eigen::VectorXd to_monomial(
			const Eigen::Index n,
			const Eigen::VectorXd& coefs)
		{
			assert((n >= 1) && "to_monomial: n must be a positive integer");

			Eigen::VectorXd coefsout = Eigen::VectorXd::Zero(n);
			Eigen::VectorXd z = Eigen::VectorXd::Zero(n);

			double v = 0.0, w = 0.0;
			coefsout(0) = coefs(n - 1);
			for (Eigen::Index j = n - 2; j >= 0; --j)
			{
				w = coefsout(0);
				coefsout(0) = coefs(j) - v * z(0);
				z(0) = w;
				for (Eigen::Index i = 1; i <= n - j - 1; ++i)
				{
					w = coefsout(i);
					coefsout(i) = (2.0 * j + 1.0) * z(i - 1) / (j + 1.0) - v * z(i);
					z(i) = w;
				}
				v = j / (j + 1.0);
			}
			return coefsout;
		}

		Eigen::VectorXd from_monomial(
			const Eigen::Index n,
			const Eigen::VectorXd& coefs)
		{
			assert((n >= 1) && "from_monomial: n must be a positive integer");

			Eigen::VectorXd coefsout = Eigen::VectorXd::Zero(n);
			Eigen::VectorXd q = Eigen::VectorXd::Zero(n);

			coefsout(0) = coefs(n - 2);
			coefsout(1) = coefs(n - 1);
			for (Eigen::Index k = 1; k <= n - 2; ++k)
			{
				q(0) = coefsout(0);
				coefsout(0) = coefs(n - k - 2) + coefsout(1) / 3.0;
				if (k >= 2)
				{
					for (Eigen::Index j = 1; j <= k - 1; ++j)
					{
						q(j) = coefsout(j);
						coefsout(j) = (j + 1.0) * coefsout(j + 1) / (2.0 * j + 3.0) + j * q(j - 1) / (2.0 * j - 1.0);
					}
				}
				q(k) = coefsout(k);
				coefsout(k) = k * q(k - 1) / (2.0 * k - 1.0);
				coefsout(k + 1) = (k + 1.0) * q(k) / (2.0 * k + 1.0);
			}
			return coefsout;
		}
	}
}
