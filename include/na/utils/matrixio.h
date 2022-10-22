#ifndef NA_UTILS_MATRIXIO_H
#define NA_UTILS_MATRIXIO_H
#include <fstream>
#include <sstream>
#include <string>
#include <typeinfo>
#include "Eigen/Core"

namespace na
{
	template <typename Derived>
	bool read_mtx(
		const std::string file_name,
		const bool match_type,
		Eigen::MatrixBase<Derived>& M)
	{
		std::ifstream file(file_name);
		if (file.fail())
		{
			return false;
		}
		std::string line;
		{
			std::getline(file, line);
			std::string beginner, object, format, field, symmetry;
			std::stringstream(line) >> beginner >> object >> format >> field >> symmetry;
			// std::string::compare() returns 0 if the two strings match exactly.
			if ((beginner.compare("%MatrixMarket") * beginner.compare("%%MatrixMarket")) != 0)
			{
				return false;
			}
			// TODO: another legal value is "vector".
			if (object.compare("matrix") != 0)
			{
				return false;
			}
			// TODO: another legal value is "coordinate", but would be more suitable to use Eigen::SparseMatrix for storage.
			if (format.compare("array") != 0)
			{
				return false;
			}
			// TODO: other legal values include "complex" and "pattern".
			if (match_type)
			{
				if (field.compare("real") == 0)
				{
					if (typeid(Derived::Scalar) != typeid(double))
					{
						return false;
					}
				}
				else if (field.compare("double") == 0)
				{
					if (typeid(Derived::Scalar) != typeid(double))
					{
						return false;
					}
				}
				else if (field.compare("integer") == 0)
				{
					if (typeid(Derived::Scalar) != typeid(int))
					{
						return false;
					}
				}
				else
				{
					return false;
				}
			}
			// TODO: other legal values include "symmetric", "skew-symmetric" and "hermitian".
			if (symmetry.compare("general") != 0)
			{
				return false;
			}
		}
		// Skip all the comment lines.
		while (std::getline(file, line))
		{
			if (line.rfind("%", 0) != 0)
			{
				break;
			}
		}
		{
			int rows, cols;
			std::stringstream(line) >> rows >> cols;
			if ((rows <= 0) || (cols <= 0))
			{
				return false;
			}
			if ((Derived::Options & Eigen::RowMajor) == 1)
			{
				M.derived().resize(cols, rows);
			}
			else
			{
				M.derived().resize(rows, cols);
			}
		}
		{
			Derived::Scalar* ptr = M.derived().data();
			while (std::getline(file, line))
			{
				std::stringstream stream(line);
				stream >> *ptr;
				if (stream.fail())
				{
					return false;
				}
				ptr += 1;
			}
			if (ptr < M.derived().data() + M.size())
			{
				return false;
			}
		}
		if ((Derived::Options & Eigen::RowMajor) == 1)
		{
			M.derived().transposeInPlace();
		}
		return true;
	}
}

#endif // !NA_UTILS_MATRIXIO_H
