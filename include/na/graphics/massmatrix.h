#ifndef NA_GRAPHICS_MASSMATRIX_H
#define NA_GRAPHICS_MASSMATIRX_H
#include <cmath>
#include "Eigen/Core"
#include "Eigen/Sparse"

namespace na
{
	namespace graphics
	{
		template <typename DerivedV, typename DerivedF, typename Scalar, int Options, typename StorageIndex>
		void massmatrix_tri_consistent(
			const Eigen::MatrixBase<DerivedV>& V,
			const Eigen::MatrixBase<DerivedF>& F,
			Eigen::SparseMatrix<Scalar, Options, StorageIndex>& M)
		{
			static_assert(DerivedV::ColsAtCompileTime == 3, "DerivedV must have compatible size");
			static_assert(DerivedF::ColsAtCompileTime == 3, "DerivedF must have compatible size");
			static_assert(std::is_same_v<typename DerivedV::Scalar, Scalar>, "V and M must have the same numeric type");
			const Eigen::Index Vrows = V.rows();
			const Eigen::Index Frows = F.rows();
			M.resize(Vrows, Vrows);
			std::vector<Eigen::Triplet<Scalar, StorageIndex>> triplets;
			triplets.reserve(9 * Frows);
			for (Eigen::Index i = 0; i < Frows; ++i)
			{
				Eigen::Vector<typename DerivedF::Scalar, 3> tri = F.row(i);
				Scalar area = static_cast<Scalar>(1.0 / 24.0) * (V.row(tri(1)) - V.row(tri(0))).cross(V.row(tri(2)) - V.row(tri(0))).norm();
				triplets.emplace_back(tri(0), tri(0), area);
				triplets.emplace_back(tri(0), tri(1), area);
				triplets.emplace_back(tri(0), tri(2), area);
				triplets.emplace_back(tri(1), tri(0), area);
				triplets.emplace_back(tri(1), tri(1), area);
				triplets.emplace_back(tri(1), tri(2), area);
				triplets.emplace_back(tri(2), tri(0), area);
				triplets.emplace_back(tri(2), tri(1), area);
				triplets.emplace_back(tri(2), tri(2), area);
			}
			M.setFromTriplets(triplets.begin(), triplets.end());
		}

		template <typename DerivedV, typename DerivedF, typename Scalar, int Options, typename StorageIndex>
		void massmatrix_tri_lumped(
			const Eigen::MatrixBase<DerivedV>& V,
			const Eigen::MatrixBase<DerivedF>& F,
			Eigen::SparseMatrix<Scalar, Options, StorageIndex>& M)
		{
			static_assert(DerivedV::ColsAtCompileTime == 3, "DerivedV must have compatible size");
			static_assert(DerivedF::ColsAtCompileTime == 3, "DerivedF must have compatible size");
			static_assert(std::is_same_v<typename DerivedV::Scalar, Scalar>, "V and M must have the same numeric type");
			const Eigen::Index Vrows = V.rows();
			const Eigen::Index Frows = F.rows();
			M.resize(Vrows, Vrows);
			M.data().reserve(Vrows);
			Eigen::Map<Eigen::Vector<Scalar, Eigen::Dynamic>, Eigen::Aligned16>(M.valuePtr(), Vrows).setZero();
			Eigen::Map<Eigen::Vector<StorageIndex, Eigen::Dynamic>, Eigen::Aligned16>(M.outerIndexPtr(), Vrows + 1).setLinSpaced(0, Vrows);
			Eigen::Map<Eigen::Vector<StorageIndex, Eigen::Dynamic>, Eigen::Aligned16>(M.innerIndexPtr(), Vrows).setLinSpaced(0, Vrows - 1);
			for (Eigen::Index i = 0; i < Frows; ++i)
			{
				Eigen::Vector<typename DerivedF::Scalar, 3> tri = F.row(i);
				Scalar area = static_cast<Scalar>(1.0 / 6.0) * (V.row(tri(1)) - V.row(tri(0))).cross(V.row(tri(2)) - V.row(tri(0))).norm();
				M.valuePtr()[tri(0)] += area;
				M.valuePtr()[tri(1)] += area;
				M.valuePtr()[tri(2)] += area;
			}
			M.data().resize(Vrows);
		}

		template <typename DerivedV, typename DerivedT, typename Scalar, int Options, typename StorageIndex>
		void massmatrix_tet(
			const Eigen::MatrixBase<DerivedV>& V,
			const Eigen::MatrixBase<DerivedT>& T,
			Eigen::SparseMatrix<Scalar, Options, StorageIndex>& M)
		{
			static_assert(DerivedV::ColsAtCompileTime == 3, "DerivedV must have compatible size");
			static_assert(DerivedT::ColsAtCompileTime == 4, "DerivedT must have compatible size");
			static_assert(std::is_same_v<typename DerivedV::Scalar, Scalar>, "V and M must have the same numeric type");
			const Eigen::Index Vrows = V.rows();
			const Eigen::Index Trows = T.rows();
			M.resize(Vrows, Vrows);
			M.data().reserve(Vrows);
			Eigen::Map<Eigen::Vector<Scalar, Eigen::Dynamic>, Eigen::Aligned16>(M.valuePtr(), Vrows).setZero();
			Eigen::Map<Eigen::Vector<StorageIndex, Eigen::Dynamic>, Eigen::Aligned16>(M.outerIndexPtr(), Vrows + 1).setLinSpaced(0, Vrows);
			Eigen::Map<Eigen::Vector<StorageIndex, Eigen::Dynamic>, Eigen::Aligned16>(M.innerIndexPtr(), Vrows).setLinSpaced(0, Vrows - 1);
			for (Eigen::Index i = 0; i < Trows; ++i)
			{
				Eigen::Vector<typename DerivedT::Scalar, 4> tet = T.row(i);
				Scalar vol = static_cast<Scalar>(1.0 / 24.0) * std::abs((V.row(tet(0)) - V.row(tet(3))).dot((V.row(tet(1)) - V.row(tet(3))).cross(V.row(tet(2)) - V.row(tet(3)))));
				M.valuePtr()[tet(0)] += vol;
				M.valuePtr()[tet(1)] += vol;
				M.valuePtr()[tet(2)] += vol;
				M.valuePtr()[tet(3)] += vol;
			}
			M.data().resize(Vrows);
		}
	}
}

#endif // !NA_GRAPHICS_MASSMATRIX_H
