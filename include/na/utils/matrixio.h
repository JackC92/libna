#ifndef NA_UTILS_MATRIXIO_H
#define NA_UTILS_MATRIXIO_H
#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "na/macros.h"
#include "na/type_traits/is_complex.h"

namespace na
{
	namespace internal
	{
		template <typename T, bool coordinate>
		constexpr const char* mtx_data_format = nullptr;

		EXPLICIT_SPEC constexpr const char* mtx_data_format<int, false> = "%i";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<float, false> = "%g";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<double, false> = "%lg";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<long double, false> = "%Lg";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<std::complex<float>, false> = "%g %g";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<std::complex<double>, false> = "%lg %lg";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<std::complex<long double>, false> = "%Lg %Lg";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<float, true> = "%i %i %g";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<double, true> = "%i %i %lg";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<long double, true> = "%i %i %Lg";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<std::complex<float>, true> = "%i %i %g %g";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<std::complex<double>, true> = "%i %i %lg %lg";
		EXPLICIT_SPEC constexpr const char* mtx_data_format<std::complex<long double>, true> = "%i %i %Lg %Lg";
	}

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
		if (!(((format.compare("array") == 0) && std::is_same_v<Derived, Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime, Derived::Options, Derived::MaxRowsAtCompileTime, Derived::MaxColsAtCompileTime>>) || ((format.compare("coordinate") == 0) && std::is_same_v<Derived, Eigen::SparseMatrix<typename Derived::Scalar, Derived::Options, typename Derived::StorageIndex>>)))
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
		if constexpr (std::is_same_v<Derived, Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime, Derived::Options, Derived::MaxRowsAtCompileTime, Derived::MaxColsAtCompileTime>>)
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

			typename Derived::RealScalar* ptr = (typename Derived::RealScalar*)M.derived().data();
			typename Derived::RealScalar* end = (typename Derived::RealScalar*)(M.derived().data() + M.size());
			while (std::getline(file, line) && (ptr < end))
			{
				if constexpr (na::is_complex_v<Derived::Scalar>)
				{
					if (std::sscanf(line.c_str(), na::internal::mtx_data_format<Derived::Scalar, false>, ptr, ptr + 1) != 2)
					{
						return false;
					}
					ptr += 2;
				}
				else
				{
					if (std::sscanf(line.c_str(), na::internal::mtx_data_format<Derived::Scalar, false>, ptr) != 1)
					{
						return false;
					}
					ptr += 1;
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

			std::vector<Eigen::Triplet<typename Derived::Scalar>> triplets;
			triplets.reserve(nnz);
			while (std::getline(file, line))
			{
				int row, col;
				typename Derived::RealScalar real, imag;
				if constexpr (na::is_complex_v<Derived::Scalar>)
				{
					if (std::sscanf(line.c_str(), na::internal::mtx_data_format<Derived::Scalar, true>, &row, &col, &real, &imag) != 4)
					{
						return false;
					}
					triplets.emplace_back(row - 1, col - 1, Derived::Scalar(real, imag));
				}
				else
				{
					if (std::sscanf(line.c_str(), na::internal::mtx_data_format<Derived::Scalar, true>, &row, &col, &real) != 3)
					{
						return false;
					}
					triplets.emplace_back(row - 1, col - 1, real);
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
