#ifndef NA_GRAPHICS_COT_LAPLACIAN_H
#define NA_GRAPHICS_COT_LAPLACIAN_H
#include <cmath>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Sparse"

namespace na
{
	namespace graphics
	{
		template <
			typename DerivedV,
			typename DerivedT,
			typename Scalar,
			int Options,
			typename StorageIndex,
			std::enable_if_t<DerivedT::ColsAtCompileTime == 4, bool> = true>
		void cot_laplacian(
			const Eigen::MatrixBase<DerivedV>& V,
			const Eigen::MatrixBase<DerivedT>& T,
			Eigen::SparseMatrix<Scalar, Options, StorageIndex>& L)
		{
			static_assert(DerivedV::ColsAtCompileTime == 3, "cot_laplacian: DerivedV must have compatible size");
			static_assert(std::is_same_v<typename DerivedV::Scalar, Scalar>, "cot_laplacian: V and L must have the same numeric type");
			const Eigen::Index Vrows = V.rows();
			const Eigen::Index Trows = T.rows();

			std::vector<Eigen::Triplet<Scalar, typename DerivedT::Scalar>> triplets;
			triplets.reserve(Trows * 16);

			Eigen::Matrix<Scalar, 3, 4> G;
			Eigen::Matrix<Scalar, 4, 4> GTG;
			for (Eigen::Index i = 0; i < Trows; ++i)
			{
				Eigen::Vector4<typename DerivedT::Scalar> tet = T.row(i);
				Eigen::Vector3<Scalar> a = V.row(tet(0));
				Eigen::Vector3<Scalar> b = V.row(tet(1));
				Eigen::Vector3<Scalar> c = V.row(tet(2));
				Eigen::Vector3<Scalar> d = V.row(tet(3));
				Eigen::Vector3<Scalar> na = (d - b).cross(c - b).normalized();
				Eigen::Vector3<Scalar> nb = (d - c).cross(a - c).normalized();
				Eigen::Vector3<Scalar> nc = (d - a).cross(b - a).normalized();
				Eigen::Vector3<Scalar> nd = (b - a).cross(c - a).normalized();

				G.col(0) = na / na.dot(a - b);
				G.col(1) = nb / nb.dot(b - c);
				G.col(2) = nc / nc.dot(c - a);
				G.col(3) = nd / nd.dot(d - a);
				GTG = (static_cast<Scalar>(1.0 / 6.0) * std::abs((V.row(tet(0)) - V.row(tet(3))).dot((V.row(tet(1)) - V.row(tet(3))).cross(V.row(tet(2)) - V.row(tet(3)))))) * (G.transpose() * G);

				for (Eigen::Index k = 0; k < 4; ++k)
				{
					triplets.emplace_back(tet(0), tet(k), GTG(0, k));
					triplets.emplace_back(tet(1), tet(k), GTG(1, k));
					triplets.emplace_back(tet(2), tet(k), GTG(2, k));
					triplets.emplace_back(tet(3), tet(k), GTG(3, k));
				}
			}

			L.resize(Vrows, Vrows);
			L.setFromTriplets(triplets.begin(), triplets.end());
		}
	}
}

#endif // !NA_GRAPHICS_COT_LAPLACIAN_H
