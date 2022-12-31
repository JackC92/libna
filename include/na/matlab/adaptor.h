#if !defined(NA_MATLAB_ADAPTOR_H) && defined(NA_USE_MATLAB)
#define NA_MATLAB_ADAPTOR_H
#include <algorithm>
#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "na/type_traits/matlab_traits.h"

namespace na
{
	namespace matlab
	{
		template <typename Derived>
		::matlab::data::TypedArray<typename Derived::Scalar> to_matlab(const Eigen::DenseBase<Derived>& mat)
		{
			typedef Eigen::Reshaped<const Derived, Derived::SizeAtCompileTime, 1> ReshapedType;
			ReshapedType vec = mat.reshaped();
			return ::matlab::data::ArrayFactory().createArray<typename ReshapedType::const_iterator, typename Derived::Scalar>({ static_cast<std::size_t>(mat.rows()), static_cast<std::size_t>(mat.cols()) }, vec.begin(), vec.end(), ::matlab::data::InputLayout::COLUMN_MAJOR);
		}

		template <typename Derived>
		::matlab::data::TypedArray<typename Derived::Scalar> to_matlab(const Eigen::PlainObjectBase<Derived>& mat)
		{
			return ::matlab::data::ArrayFactory().createArray<const typename Derived::Scalar*, typename Derived::Scalar>({ static_cast<std::size_t>(mat.rows()), static_cast<std::size_t>(mat.cols()) }, mat.data(), mat.data() + mat.size(), static_cast<::matlab::data::InputLayout>(Derived::IsRowMajor));
		}

		template <typename Scalar, int Options, typename StorageIndex>
		::matlab::data::SparseArray<typename na::matlab::sparse_buffer_traits<Scalar>::data_type> to_matlab(const Eigen::SparseMatrix<Scalar, Options, StorageIndex>& mat)
		{
			typedef typename na::matlab::sparse_buffer_traits<Scalar>::data_type ArrayScalar;

			const std::size_t nnz = mat.nonZeros();
			::matlab::data::ArrayFactory factory;
			::matlab::data::buffer_ptr_t<ArrayScalar> data = factory.createBuffer<ArrayScalar>(nnz);
			::matlab::data::buffer_ptr_t<std::size_t> outer = factory.createBuffer<std::size_t>(nnz);
			::matlab::data::buffer_ptr_t<std::size_t> inner = factory.createBuffer<std::size_t>(nnz);
			if (!mat.isCompressed())
			{
				ArrayScalar* data_ptr = data.get();
				std::size_t* inner_ptr = inner.get();
				std::size_t* outer_ptr = outer.get();
				for (Eigen::Index i = 0; i < mat.outerSize(); ++i)
				{
					const Eigen::Index inner_nnz = mat.innerVector(i).nonZeros();
					Eigen::Map<Eigen::VectorX<ArrayScalar>>(data_ptr, inner_nnz) = Eigen::Map<const Eigen::VectorX<Scalar>>(mat.valuePtr() + mat.outerIndexPtr()[i], inner_nnz).template cast<ArrayScalar>();
					Eigen::Map<Eigen::VectorX<std::size_t>>(inner_ptr, inner_nnz) = Eigen::Map<const Eigen::VectorX<StorageIndex>>(mat.innerIndexPtr() + mat.outerIndexPtr()[i], inner_nnz).template cast<std::size_t>();
					Eigen::Map<Eigen::VectorX<std::size_t>>(outer_ptr, inner_nnz).setConstant(static_cast<std::size_t>(i));
					data_ptr += inner_nnz;
					inner_ptr += inner_nnz;
					outer_ptr += inner_nnz;
				}

			}
			else
			{
				Eigen::Map<Eigen::VectorX<ArrayScalar>>(data.get(), nnz) = Eigen::Map<const Eigen::VectorX<Scalar>>(mat.valuePtr(), nnz).template cast<ArrayScalar>();
				Eigen::Map<Eigen::VectorX<std::size_t>>(inner.get(), nnz) = Eigen::Map<const Eigen::VectorX<StorageIndex>>(mat.innerIndexPtr(), nnz).template cast<std::size_t>();
				for (Eigen::Index i = 0; i < mat.outerSize(); ++i)
				{
					Eigen::Map<Eigen::VectorX<std::size_t>>(outer.get() + mat.outerIndexPtr()[i], mat.innerVector(i).nonZeros()).setConstant(static_cast<std::size_t>(i));
				}
			}
			if constexpr (Eigen::SparseMatrix<Scalar, Options, StorageIndex>::IsRowMajor)
			{
				return factory.createSparseArray<ArrayScalar>({ static_cast<std::size_t>(mat.rows()), static_cast<std::size_t>(mat.cols()) }, nnz, std::move(data), std::move(outer), std::move(inner));
			}
			else
			{
				return factory.createSparseArray<ArrayScalar>({ static_cast<std::size_t>(mat.rows()), static_cast<std::size_t>(mat.cols()) }, nnz, std::move(data), std::move(inner), std::move(outer));
			}
		}

		template <typename Derived>
		void to_eigen(
			const ::matlab::data::TypedArrayRef<typename Derived::Scalar>& arr,
			Eigen::PlainObjectBase<Derived>& mat)
		{
			mat.resize(static_cast<Eigen::Index>(arr.getDimensions()[0]), static_cast<Eigen::Index>(arr.getDimensions()[1]));
			for (Eigen::Index j = 0; j < mat.cols(); ++j)
			{
				for (Eigen::Index i = 0; i < mat.rows(); ++i)
				{
					mat.coeffRef(i, j) = arr[i][j];
				}
			}
		}

		template <typename Derived>
		void to_eigen(
			const ::matlab::data::TypedArrayRef<typename Derived::Scalar>& arr,
			Eigen::DiagonalBase<Derived>& mat)
		{
			mat.diagonal().resize(static_cast<Eigen::Index>(std::min<std::size_t>(arr.getDimensions()[0], arr.getDimensions()[1])));
			for (Eigen::Index i = 0; i < mat.rows(); ++i)
			{
				mat.diagonal()(i) = arr[i][i];
			}
		}
	}
}

#endif // !defined(NA_MATLAB_ADAPTOR_H) && defined(NA_USE_MATLAB)
