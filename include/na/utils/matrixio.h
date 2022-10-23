#ifndef NA_UTILS_MATRIXIO_H
#define NA_UTILS_MATRIXIO_H
#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "na/type_traits/is_complex.h"

namespace na
{
	template <typename Derived>
	bool read_mtx(
		const std::string& file_name,
		Eigen::EigenBase<Derived>& M)
	{
		std::ifstream file(file_name);
		if (file.fail())
		{
			return false;
		}
		std::string line;
		std::getline(file, line);
		std::string beginner, object, format, field, symmetry;
		std::stringstream(line) >> beginner >> object >> format >> field >> symmetry;
		// std::string::compare() returns 0 if the two strings match exactly.
		if ((beginner.compare("%MatrixMarket") * beginner.compare("%%MatrixMarket")) != 0)
		{
			return false;
		}
		if (object.compare("matrix") != 0)
		{
			return false;
		}
		// This establishes the equivalence between "array" and Eigen::Matrix, or "coordinate" and Eigen::SparseMatrix.
		if (!(((format.compare("array") == 0) && std::is_same_v<Derived, Eigen::Matrix<Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime, Derived::Options, Derived::MaxRowsAtCompileTime, Derived::MaxColsAtCompileTime>>) || ((format.compare("coordinate") == 0) && std::is_same_v<Derived, Eigen::SparseMatrix<Derived::Scalar, Derived::Options, Derived::StorageIndex>>)))
		{
			return false;
		}
		if ((format.compare("coordinate") == 0) ^ (field.compare("pattern") == 0))
		{
			return false;
		}
		// TODO: other legal values include "symmetric", "skew-symmetric" and "hermitian".
		if (symmetry.compare("general") != 0)
		{
			return false;
		}
		// Skip all the comment lines.
		while (std::getline(file, line))
		{
			if (line.rfind("%", 0) != 0)
			{
				break;
			}
		}
		// The while loop above should have consumed the size line, if it exists. The following code does the actual check on the size line.
		if constexpr (std::is_same_v<Derived, Eigen::Matrix<Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime, Derived::Options, Derived::MaxRowsAtCompileTime, Derived::MaxColsAtCompileTime>>)
		{
			int rows, cols;
			std::stringstream(line) >> rows >> cols;
			if ((rows <= 0) || (cols <= 0))
			{
				return false;
			}
			if ((Derived::RowsAtCompileTime != Eigen::Dynamic) && (rows != Derived::RowsAtCompileTime))
			{
				return false;
			}
			if ((Derived::ColsAtCompileTime != Eigen::Dynamic) && (cols != Derived::ColsAtCompileTime))
			{
				return false;
			}
			if ((Derived::RowsAtCompileTime == Eigen::Dynamic) && (Derived::MaxRowsAtCompileTime != Eigen::Dynamic) && (rows > Derived::MaxRowsAtCompileTime))
			{
				return false;
			}
			if ((Derived::ColsAtCompileTime == Eigen::Dynamic) && (Derived::MaxColsAtCompileTime != Eigen::Dynamic) && (cols > Derived::MaxColsAtCompileTime))
			{
				return false;
			}
			M.derived().resize(rows, cols);

			Derived::RealScalar* ptr = (Derived::RealScalar*)M.derived().data();
			Derived::RealScalar* end = (Derived::RealScalar*)(M.derived().data() + M.size());
			while (std::getline(file, line) && (ptr < end))
			{
				std::stringstream stream(line);
				if constexpr (na::is_complex_v<Derived::Scalar>)
				{
					stream >> *ptr >> *(ptr + 1);
					ptr += 2;
				}
				else
				{
					stream >> *ptr;
					ptr += 1;
				}
				if (stream.fail())
				{
					return false;
				}
			}
			if (ptr < end)
			{
				return false;
			}
			if ((Derived::Options & Eigen::RowMajor) == 1)
			{
				M.derived() = M.derived().transpose().reshaped(rows, cols).eval();
			}
		}
		else
		{
			int rows, cols, nnz;
			std::stringstream(line) >> rows >> cols >> nnz;
			if ((rows <= 0) || (cols <= 0) || (nnz <= 0))
			{
				return false;
			}
			M.derived().resize(rows, cols);
			M.derived().reserve(nnz);

			std::vector<Eigen::Triplet<Derived::Scalar>> triplets;
			triplets.reserve(nnz);
			while (std::getline(file, line))
			{
				int row, col;
				Derived::RealScalar real, imag;
				std::stringstream stream(line);
				stream >> row >> col;
				if constexpr (na::is_complex_v<Derived::Scalar>)
				{
					stream >> real >> imag;
					triplets.emplace_back(row - 1, col - 1, Derived::Scalar(real, imag));
				}
				else
				{
					stream >> real;
					triplets.emplace_back(row - 1, col - 1, real);
				}
				if (stream.fail())
				{
					return false;
				}
			}
			if (triplets.size() != nnz)
			{
				return false;
			}
			M.derived().setFromTriplets(triplets.begin(), triplets.end());
		}
		return true;
	}
}

#endif // !NA_UTILS_MATRIXIO_H
