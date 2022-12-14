#ifndef NA_IO_MTX_H
#define NA_IO_MTX_H
#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "na/macros.h"
#include "na/type_traits/is_array_or_matrix.h"
#include "na/type_traits/is_complex.h"

namespace na
{
	namespace detail
	{
		template <typename T, bool coordinate>
		constexpr const char* mtx_data_format = nullptr;

		EXPSPEC constexpr const char* mtx_data_format<int, false> = "%d";
		EXPSPEC constexpr const char* mtx_data_format<float, false> = "%g";
		EXPSPEC constexpr const char* mtx_data_format<double, false> = "%lg";
		EXPSPEC constexpr const char* mtx_data_format<long double, false> = "%Lg";
		EXPSPEC constexpr const char* mtx_data_format<std::complex<float>, false> = "%g %g";
		EXPSPEC constexpr const char* mtx_data_format<std::complex<double>, false> = "%lg %lg";
		EXPSPEC constexpr const char* mtx_data_format<std::complex<long double>, false> = "%Lg %Lg";
		EXPSPEC constexpr const char* mtx_data_format<float, true> = "%lld %lld %g";
		EXPSPEC constexpr const char* mtx_data_format<double, true> = "%lld %lld %lg";
		EXPSPEC constexpr const char* mtx_data_format<long double, true> = "%lld %lld %Lg";
		EXPSPEC constexpr const char* mtx_data_format<std::complex<float>, true> = "%lld %lld %g %g";
		EXPSPEC constexpr const char* mtx_data_format<std::complex<double>, true> = "%lld %lld %lg %lg";
		EXPSPEC constexpr const char* mtx_data_format<std::complex<long double>, true> = "%lld %lld %Lg %Lg";
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
		// This establishes the equivalence between "array" and Eigen::PlainObjectBase<Derived>, or "coordinate" and Eigen::SparseMatrix.
		if (!(((format.compare("array") == 0) && na::is_array_or_matrix_v<Derived>) || ((format.compare("coordinate") == 0) && std::is_same_v<Derived, Eigen::SparseMatrix<typename Derived::Scalar, Derived::Options, typename Derived::StorageIndex>>)))
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
		if constexpr (na::is_array_or_matrix_v<Derived>)
		{
			Eigen::Index rows, cols;
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
				if constexpr (na::is_complex_v<typename Derived::Scalar>)
				{
					if (std::sscanf(line.c_str(), na::detail::mtx_data_format<typename Derived::Scalar, false>, ptr, ptr + 1) != 2)
					{
						return false;
					}
					ptr += 2;
				}
				else
				{
					if (std::sscanf(line.c_str(), na::detail::mtx_data_format<typename Derived::Scalar, false>, ptr) != 1)
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
			Eigen::Index rows, cols, nnz;
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
				Eigen::Index row, col;
				typename Derived::RealScalar real, imag;
				if constexpr (na::is_complex_v<typename Derived::Scalar>)
				{
					if (std::sscanf(line.c_str(), na::detail::mtx_data_format<typename Derived::Scalar, true>, &row, &col, &real, &imag) != 4)
					{
						return false;
					}
					triplets.emplace_back(row - 1, col - 1, Derived::Scalar(real, imag));
				}
				else
				{
					if (std::sscanf(line.c_str(), na::detail::mtx_data_format<typename Derived::Scalar, true>, &row, &col, &real) != 3)
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

#endif // !NA_IO_MTX_H
