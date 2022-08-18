#ifndef NA_LINALG_VECTOR_H
#define NA_LINALG_VECTOR_H
#include <type_traits>
#include "Eigen/Core"
#include "na/type_traits/floating_point_scalar.h"

namespace na
{
	namespace linalg
	{
		// TODO: extend for long double, considering both 80-bit extended precision format and 128-bit quadruple precision format.
		template <typename Real>
		struct scaling_constants
		{
			static_assert(std::is_floating_point_v<Real>, "scaling_constants can only be used for floating point types");
			static Real tau_min;
			static Real tau_max;
			static Real sig_min;
			static Real sig_max;
		};

		// This function is based on the paper
		//   Walter F. Mascarenhas, "Fast and accurate normalization of vectors and quaternions"
		template <typename Scalar, int Rows, typename Real = na::floating_point_scalar<Scalar>::value_type>
		inline void scale(
			Eigen::Vector<Scalar, Rows>& x,
			Real& coeff)
		{
			Real m = x.cwiseAbs().maxCoeff();
			if (m == static_cast<Real>(0.0))
			{
				coeff = static_cast<Real>(0.0);
				return;
			}
			if (m >= scaling_constants<Real>::tau_min)
			{
				if (m <= scaling_constants<Real>::tau_max)
				{
					coeff = static_cast<Real>(1.0);
				}
				else
				{
					x *= scaling_constants<Real>::sig_max;
					coeff = static_cast<Real>(1.0) / scaling_constants<Real>::sig_max;
				}
			}
			else
			{
				x *= scaling_constants<Real>::sig_min;
				coeff = static_cast<Real>(1.0) / scaling_constants<Real>::sig_min;
			}
		}

		template <typename Scalar, int Rows, typename Real = na::floating_point_scalar<Scalar>::value_type>
		inline void scale(Eigen::Vector<Scalar, Rows>& x)
		{
			Real m = x.cwiseAbs().maxCoeff();
			if (m >= scaling_constants<Real>::tau_min)
			{
				if (m > scaling_constants<Real>::tau_max)
				{
					x *= scaling_constants<Real>::sig_max;
				}
			}
			else
			{
				x *= scaling_constants<Real>::sig_min;
			}
		}

		template <typename Scalar, int Rows, typename Real = na::floating_point_scalar<Scalar>::value_type>
		inline Real norm(const Eigen::Vector<Scalar, Rows>& x)
		{
			Real coeff;
			Eigen::Vector<Scalar, Rows> vec = x;
			scale(vec, coeff);
			return coeff * vec.norm();
		}

		template <typename Scalar, int Rows>
		inline void normalize(Eigen::Vector<Scalar, Rows>& x)
		{
			scale(x);
			x.normalize();
		}

		template <typename Scalar, int Rows>
		inline Eigen::Vector<Scalar, Rows> normalized(const Eigen::Vector<Scalar, Rows>& x)
		{
			Eigen::Vector<Scalar, Rows> vec = x;
			normalize(vec);
			return vec;
		}
	}
}

#endif // !NA_LINALG_VECTOR_H
