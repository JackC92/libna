#ifndef NA_CORE_SPARSE_H
#define NA_CORE_SPARSE_H
#include <cstring>
#include <numeric>
#include "Eigen/Core"
#include "Eigen/Sparse"

namespace na
{
	template <typename Scalar, int Options, typename StorageIndex>
	inline void repmat(
		const Eigen::SparseMatrix<Scalar, Options, StorageIndex>& mat,
		const Eigen::Index m,
		const Eigen::Index n,
		Eigen::SparseMatrix<Scalar, Options, StorageIndex>& res)
	{
		typedef Eigen::SparseMatrix<Scalar, Options, StorageIndex> MatrixType;

		const Eigen::Index nnz = mat.nonZeros();
		const Eigen::Index outer_factor = MatrixType::IsRowMajor ? m : n;
		const Eigen::Index inner_factor = MatrixType::IsRowMajor ? n : m;

		res.resize(mat.rows() * m, mat.cols() * n);
		res.data().reserve(nnz * m * n);
		Scalar* value_ptr = res.valuePtr();
		StorageIndex* inner_ptr = res.innerIndexPtr();
		if (!mat.isCompressed())
		{
			std::partial_sum<const StorageIndex*, StorageIndex*>(mat.innerNonZeroPtr(), mat.innerNonZeroPtr() + mat.outerSize() - 1, res.outerIndexPtr() + 1);
			Eigen::Map<Eigen::ArrayX<StorageIndex>>(res.outerIndexPtr(), mat.outerSize()) *= inner_factor;
		}
		else
		{
			Eigen::Map<Eigen::ArrayX<StorageIndex>>(res.outerIndexPtr(), mat.outerSize()) = inner_factor * Eigen::Map<const Eigen::ArrayX<StorageIndex>>(mat.outerIndexPtr(), mat.outerSize());
		}
		for (Eigen::Index i = 0; i < mat.outerSize(); ++i)
		{
			const Eigen::Index inner_nnz = mat.innerVector(i).nonZeros();
			for (Eigen::Index j = 0; j < inner_factor; ++j)
			{
				std::memcpy(value_ptr, mat.valuePtr() + mat.outerIndexPtr()[i], sizeof(Scalar) * inner_nnz);
				Eigen::Map<Eigen::ArrayX<StorageIndex>>(inner_ptr, inner_nnz) = Eigen::Map<const Eigen::ArrayX<StorageIndex>>(mat.innerIndexPtr() + mat.outerIndexPtr()[i], inner_nnz) + j * mat.innerSize();
				value_ptr += inner_nnz;
				inner_ptr += inner_nnz;
			}
		}
		for (Eigen::Index i = 1; i < outer_factor; ++i)
		{
			Eigen::Map<Eigen::ArrayX<StorageIndex>>(res.outerIndexPtr() + i * mat.outerSize(), mat.outerSize()) = Eigen::Map<const Eigen::ArrayX<StorageIndex>>(res.outerIndexPtr(), mat.outerSize()) + i * inner_factor * nnz;
			std::memcpy(value_ptr, res.valuePtr(), sizeof(Scalar) * inner_factor * nnz);
			std::memcpy(inner_ptr, res.innerIndexPtr(), sizeof(StorageIndex) * inner_factor * nnz);
			value_ptr += inner_factor * nnz;
			inner_ptr += inner_factor * nnz;
		}
		res.outerIndexPtr()[res.outerSize()] = nnz * m * n;
		res.data().resize(nnz * m * n);
	}

	template <typename Scalar, int Options, typename StorageIndex>
	inline void repdiag(
		const Eigen::SparseMatrix<Scalar, Options, StorageIndex>& mat,
		const Eigen::Index m,
		Eigen::SparseMatrix<Scalar, Options, StorageIndex>& res)
	{
		const Eigen::Index rows = mat.rows();
		const Eigen::Index cols = mat.cols();
		const Eigen::Index nnz = mat.nonZeros();

		res.resize(rows * m, cols * m);
		res.data().reserve(nnz * m);
		if (!mat.isCompressed())
		{
			std::partial_sum<const StorageIndex*, StorageIndex*>(mat.innerNonZeroPtr(), mat.innerNonZeroPtr() + mat.outerSize() - 1, res.outerIndexPtr() + 1);
			for (Eigen::Index i = 0; i < mat.outerSize(); ++i)
			{
				std::memcpy(res.valuePtr() + count, mat.valuePtr() + mat.outerIndexPtr()[i], sizeof(Scalar) * mat.innerNonZeroPtr()[i]);
				std::memcpy(res.innerIndexPtr() + count, mat.innerIndexPtr() + mat.outerIndexPtr()[i], sizeof(StorageIndex) * mat.innerNonZeroPtr()[i]);
			}
		}
		else
		{
			std::memcpy(res.valuePtr(), mat.valuePtr(), sizeof(Scalar) * nnz);
			std::memcpy(res.innerIndexPtr(), mat.innerIndexPtr(), sizeof(StorageIndex) * nnz);
			std::memcpy(res.outerIndexPtr(), mat.outerIndexPtr(), sizeof(StorageIndex) * mat.outerSize());
		}
		for (Eigen::Index i = 1; i < m; ++i)
		{
			std::memcpy(res.valuePtr() + nnz * i, res.valuePtr(), sizeof(Scalar) * nnz);
			std::memcpy(res.innerIndexPtr() + nnz * i, res.innerIndexPtr(), sizeof(StorageIndex) * nnz);
			std::memcpy(res.outerIndexPtr() + mat.outerSize() * i, res.outerIndexPtr(), sizeof(StorageIndex) * mat.outerSize());
			Eigen::Map<Eigen::ArrayX<StorageIndex>>(res.outerIndexPtr() + mat.outerSize() * i, mat.outerSize()) += nnz * i;
			Eigen::Map<Eigen::ArrayX<StorageIndex>>(res.innerIndexPtr() + nnz * i, nnz) += mat.innerSize() * i;
		}
		res.outerIndexPtr()[res.outerSize()] = nnz * m;
		res.data().resize(nnz * m);
	}
}

#endif // !NA_CORE_SPARSE_H
