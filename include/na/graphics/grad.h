#ifndef NA_GRAPHICS_GRAD_H
#define NA_GRAPHICS_GRAD_H
#include <type_traits>
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/Sparse"

namespace na
{
	namespace graphics
	{
		// Compute the per-tetrahedron gradient of a piecewise-linear scalar function defined on the vertices.
		// Note that the gradient is constant over each tetrahedron.
		template <
			typename DerivedV,
			typename DerivedT,
			typename Scalar,
			typename StorageIndex,
			std::enable_if_t<DerivedT::ColsAtCompileTime == 4, bool> = true>
		void grad(
			const Eigen::MatrixBase<DerivedV>& V,
			const Eigen::MatrixBase<DerivedT>& T,
			Eigen::SparseMatrix<Scalar, Eigen::RowMajor, StorageIndex>& G)
		{
			static_assert(DerivedV::ColsAtCompileTime == 3, "grad: DerivedV must have compatible size");
			static_assert(std::is_same_v<typename DerivedV::Scalar, Scalar>, "grad: V and G must have the same numeric type");
			const Eigen::Index Vrows = V.rows();
			const Eigen::Index Trows = T.rows();
			// This updates the size of the matrix and reallocates memory for the outer-index array.
			G.resize(3 * Trows, Vrows);
			// This reallocates sufficient memory for the value and inner-index array, but does not update the size tracked by the CompressedStorage object.
			G.data().reserve(4 * G.rows());
			Eigen::Map<Eigen::VectorX<StorageIndex>, Eigen::Aligned16>(G.outerIndexPtr(), 3 * Trows + 1).setLinSpaced(0, 12 * Trows);
			Eigen::Vector4<StorageIndex> order(0, 1, 2, 3), argorder(0, 1, 2, 3);
			for (Eigen::Index i = 0; i < Trows; ++i)
			{
				Eigen::Vector4<typename DerivedT::Scalar> tet = T.row(i);
				std::sort<StorageIndex*>(order.data(), order.data() + 4, [&tet](StorageIndex k, StorageIndex l) { return tet(k) < tet(l); });
				std::sort<StorageIndex*>(argorder.data(), argorder.data() + 4, [&order](StorageIndex k, StorageIndex l) { return order(k) < order(l); });
				Eigen::Vector3<Scalar> a = V.row(tet(0));
				Eigen::Vector3<Scalar> b = V.row(tet(1));
				Eigen::Vector3<Scalar> c = V.row(tet(2));
				Eigen::Vector3<Scalar> d = V.row(tet(3));
				Eigen::Vector3<Scalar> na = (d - b).cross(c - b).normalized();
				Eigen::Vector3<Scalar> nb = (d - c).cross(a - c).normalized();
				Eigen::Vector3<Scalar> nc = (d - a).cross(b - a).normalized();
				Eigen::Vector3<Scalar> nd = (b - a).cross(c - a).normalized();
				Eigen::Map<Eigen::Vector3<Scalar>, Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>>(G.valuePtr() + 4 * i + argorder(0), Eigen::InnerStride<Eigen::Dynamic>(4 * Trows)) = na / na.dot(a - b);
				Eigen::Map<Eigen::Vector3<Scalar>, Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>>(G.valuePtr() + 4 * i + argorder(1), Eigen::InnerStride<Eigen::Dynamic>(4 * Trows)) = nb / nb.dot(b - c);
				Eigen::Map<Eigen::Vector3<Scalar>, Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>>(G.valuePtr() + 4 * i + argorder(2), Eigen::InnerStride<Eigen::Dynamic>(4 * Trows)) = nc / nc.dot(c - a);
				Eigen::Map<Eigen::Vector3<Scalar>, Eigen::Unaligned, Eigen::InnerStride<Eigen::Dynamic>>(G.valuePtr() + 4 * i + argorder(3), Eigen::InnerStride<Eigen::Dynamic>(4 * Trows)) = nd / nd.dot(d - a);
				Eigen::Map<Eigen::Matrix<StorageIndex, 4, 3>, Eigen::Aligned16, Eigen::OuterStride<Eigen::Dynamic>>(G.innerIndexPtr() + 4 * i, Eigen::OuterStride<Eigen::Dynamic>(4 * Trows)) << tet(order), tet(order), tet(order);
			}
			G.data().resize(4 * G.rows());
		}

		template <
			typename DerivedV,
			typename DerivedT,
			typename Scalar,
			typename StorageIndex,
			std::enable_if_t<DerivedT::ColsAtCompileTime == 4, bool> = true>
		void grad(
			const Eigen::MatrixBase<DerivedV>& V,
			const Eigen::MatrixBase<DerivedT>& T,
			Eigen::SparseMatrix<Scalar, Eigen::ColMajor, StorageIndex>& G)
		{
			Eigen::SparseMatrix<Scalar, Eigen::RowMajor, StorageIndex> trG;
			grad(V, T, trG);
			G = trG;
		}
	}
}

#endif // !NA_GRAPHICS_GRAD_H
