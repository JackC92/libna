#ifndef NA_LINALG_SPARSE_H
#define NA_LINALG_SPARSE_H
#include <cstring>
#include "Eigen/Core"
#include "Eigen/Sparse"

namespace na
{
	namespace linalg
	{
		template <typename Scalar, int Options, typename StorageIndex>
		void repdiag(
			const Eigen::SparseMatrix<Scalar, Options, StorageIndex>& A,
			const Eigen::Index d,
			Eigen::SparseMatrix<Scalar, Options, StorageIndex>& B)
		{
			const Eigen::Index rows = A.rows();
			const Eigen::Index cols = A.cols();
			const Eigen::Index nnz = A.nonZeros();

			B.resize(rows * d, cols * d);
			B.data().reserve(nnz * d);
			if (!A.isCompressed())
			{
				StorageIndex count = 0;
				for (Eigen::Index j = 0; j < A.outerSize(); ++j)
				{
					B.outerIndexPtr()[j] = count;
					std::memcpy(B.innerIndexPtr() + count, A.innerIndexPtr() + A.outerIndexPtr()[j], sizeof(StorageIndex) * A.innerNonZeroPtr()[j]);
					std::memcpy(B.valuePtr() + count, A.valuePtr() + A.outerIndexPtr()[j], sizeof(Scalar) * A.innerNonZeroPtr()[j]);
					count += A.innerNonZeroPtr()[j];
				}
			}
			else
			{
				std::memcpy(B.valuePtr(), A.valuePtr(), sizeof(Scalar) * nnz);
				std::memcpy(B.innerIndexPtr(), A.innerIndexPtr(), sizeof(StorageIndex) * nnz);
				std::memcpy(B.outerIndexPtr(), A.outerIndexPtr(), sizeof(StorageIndex) * A.outerSize());
			}
			for (Eigen::Index i = 1; i < d; ++i)
			{
				std::memcpy(B.valuePtr() + nnz * i, A.valuePtr(), sizeof(Scalar) * nnz);
				std::memcpy(B.innerIndexPtr() + nnz * i, A.innerIndexPtr(), sizeof(StorageIndex) * nnz);
				std::memcpy(B.outerIndexPtr() + A.outerSize() * i, A.outerIndexPtr(), sizeof(StorageIndex) * A.outerSize());
				Eigen::Map<Eigen::Array<StorageIndex, Eigen::Dynamic, 1>>(B.outerIndexPtr() + A.outerSize() * i, A.outerSize()) += nnz * i;
				Eigen::Map<Eigen::Array<StorageIndex, Eigen::Dynamic, 1>>(B.innerIndexPtr() + nnz * i, nnz) += A.innerSize() * i;
			}
			B.outerIndexPtr()[B.outerSize()] = nnz * d;
			B.data().resize(nnz * d);
		}
	}
}

#endif // !NA_LINALG_SPARSE_H
