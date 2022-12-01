#ifndef NA_LINALG_VECTOR_H
#define NA_LINALG_VECTOR_H
#include <algorithm>
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
			static constexpr Real tau_min{ 0.0 };
			static constexpr Real tau_max{ 0.0 };
			static constexpr Real sig_min{ 0.0 };
			static constexpr Real sig_max{ 0.0 };
		};

		template <>
		struct scaling_constants<float>
		{
			static constexpr float tau_min{ 1.7763568394002504646778106689453e-15f };
			static constexpr float tau_max{ 4.6116860184273879040000000000000e+18f };
			static constexpr float sig_min{ 1.2676506002282294014967032053760e+30f };
			static constexpr float sig_max{ 1.3552527156068805425093160010874e-20f };
		};

		template <>
		struct scaling_constants<double>
		{
			static constexpr double tau_min{ 8.0083323807324036977183408605459e-146 };
			static constexpr double tau_max{ 3.3519519824856492748935062495515e+153 };
			static constexpr double sig_min{ 1.6209045190941378744189093217544e+178 };
			static constexpr double sig_max{ 1.8645851828000516858227413288657e-155 };
		};

		// This function is based on the paper
		//   Walter F. Mascarenhas, "Fast and accurate normalization of vectors and quaternions"
		template <typename Scalar, int Rows, typename Real = typename na::floating_point_scalar<Scalar>::value_type>
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

		template <typename Scalar, int Rows, typename Real = typename na::floating_point_scalar<Scalar>::value_type>
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

		template <typename Scalar, int Rows, typename Real = typename na::floating_point_scalar<Scalar>::value_type>
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
		
		template <typename Scalar>
		inline void deduplicate(
			Eigen::Vector<Scalar, Eigen::Dynamic>& x,
			const bool sort = false)
		{
			if (x.size() == 0)
			{
				return;
			}
			if (sort)
			{
				std::sort<Scalar*>(x.data(), x.data() + x.size());
			}
			// This does not guarantee uniqueness of the entries unless the vector is sorted.
			x.conservativeResize(std::unique<Scalar*>(x.data(), x.data() + x.size()) - x.data());
		}

		template <typename Derived, typename Scalar = typename Derived::Scalar>
		inline Eigen::Vector<Scalar, Eigen::Dynamic> repelem(
			const Eigen::MatrixBase<Derived>& x,
			const Eigen::Index n)
		{
			static_assert((Derived::RowsAtCompileTime == 1) || (Derived::ColsAtCompileTime == 1), "Derived must be of the shape of a vector");
			// The first reshaped() is equivalent to Matlab's x(:), which results in a column vector.
			return x.reshaped().transpose().replicate(n, 1).reshaped();
		}
	}
}

#endif // !NA_LINALG_VECTOR_H
