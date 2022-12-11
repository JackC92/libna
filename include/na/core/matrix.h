#ifndef NA_CORE_MATRIX_H
#define NA_CORE_MATRIX_H
#include <cassert>
#include <type_traits>
#include "Eigen/Core"
#ifdef NA_USE_MKL
#include "mkl_vsl.h"
#endif // NA_USE_MKL
#include "na/mkl/distributions.h"
#include "na/type_traits/is_complex.h"

namespace na
{
	template <typename Derived>
	inline void repelem(
		const Eigen::MatrixBase<Derived>& mat,
		const Eigen::Index m,
		const Eigen::Index n,
		Eigen::MatrixX<typename Derived::Scalar>& res)
	{
		res.resize(mat.rows() * m, mat.cols() * n);
		for (Eigen::Index i = 0; i < mat.rows(); ++i)
		{
			for (Eigen::Index j = 0; j < mat.cols(); ++j)
			{
				res.block(i * m, j * n, m, n).setConstant(mat(i, j));
			}
		}
	}

	template <typename Derived, typename Index>
	inline void repelem(
		const Eigen::MatrixBase<Derived>& mat,
		const Eigen::VectorX<Index>& m,
		const Eigen::VectorX<Index>& n,
		Eigen::MatrixX<typename Derived::Scalar>& res)
	{
		assert(((mat.rows() == m.size()) && (mat.cols() == n.size())) && "repelem: m and n must have the same size as the number of rows and columns of mat");

		res.resize(mat.rows() * m.sum(), mat.cols() * n.sum());
		Eigen::Index i = 0;
		for (Eigen::Index k = 0; k < mat.rows(); ++k)
		{
			Eigen::Index j = 0;
			for (Eigen::Index l = 0; l < mat.cols(); ++l)
			{
				res.block(i, j, m(k), n(l)).setConstant(mat(k, l));
				j += n(l);
			}
			i += m(k);
		}
	}

#ifdef NA_USE_MKL
	/*
	* The following subroutines rely on Intel MKL, and are only applicable to int, float, and double because they are the data types supported by MKL.
	*/
	template <
		typename Derived,
		std::enable_if_t<std::is_floating_point_v<typename Derived::Scalar>, bool> = true>
	inline void set_random(
		Eigen::PlainObjectBase<Derived>& mat,
		VSLStreamStatePtr stream = nullptr,
		const typename Derived::Scalar lower = typename Derived::Scalar(-1.0),
		const typename Derived::Scalar upper = typename Derived::Scalar(+1.0))
	{
		VSLStreamStatePtr ptr = stream;
		const bool has_stream = (ptr != nullptr);
		if (!has_stream)
		{
			// Stream initialization is likely to be costly, so it should be prefered to pass an initialized stream to this function.
			vslNewStream(&ptr, VSL_BRNG_MT19937, 0);
		}
		// Note that the MKL's vRngUniform functions exclude the right bound of the interval.
		na::mkl::distribution_generator<typename Derived::Scalar, na::mkl::DIST_UNIFORM>::run(VSL_RNG_METHOD_UNIFORM_STD, ptr, mat.size(), mat.derived().data(), lower, upper);
		if (!has_stream)
		{
			vslDeleteStream(&ptr);
		}
	}

	template <
		typename Derived,
		std::enable_if_t<na::is_complex_v<typename Derived::Scalar>, bool> = true>
	inline void set_random(
		Eigen::PlainObjectBase<Derived>& mat,
		VSLStreamStatePtr stream = nullptr,
		const typename Derived::Scalar lower = typename Derived::Scalar(-1.0, -1.0),
		const typename Derived::Scalar upper = typename Derived::Scalar(+1.0, +1.0))
	{
		VSLStreamStatePtr ptr = stream;
		const bool has_stream = (ptr != nullptr);
		if (!has_stream)
		{
			vslNewStream(&ptr, VSL_BRNG_MT19937, 0);
		}
		na::mkl::distribution_generator<typename Derived::RealScalar, na::mkl::DIST_UNIFORM>::run(VSL_RNG_METHOD_UNIFORM_STD, ptr, mat.size(), (typename Derived::RealScalar*)mat.derived().data(), lower.real(), upper.real());
		na::mkl::distribution_generator<typename Derived::RealScalar, na::mkl::DIST_UNIFORM>::run(VSL_RNG_METHOD_UNIFORM_STD, ptr, mat.size(), ((typename Derived::RealScalar*)mat.derived().data()) + mat.size(), lower.imag(), upper.imag());
		Eigen::Map<Eigen::Matrix2X<typename Derived::RealScalar>>((typename Derived::RealScalar*)mat.derived().data(), 2, mat.size()) = Eigen::Map<Eigen::MatrixX2<typename Derived::RealScalar>>((typename Derived::RealScalar*)mat.derived().data(), mat.size(), 2).transpose().eval();
		if (!has_stream)
		{
			vslDeleteStream(&ptr);
		}
	}
#endif // NA_USE_MKL
}

#endif // !NA_CORE_MATRIX_H
