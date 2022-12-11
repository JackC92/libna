#ifndef NA_IO_MATLAB_H
#define NA_IO_MATLAB_H
#include <complex>
#include <fstream>
#include <string>
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "na/type_traits/is_complex.h"

namespace na
{
	template <typename Scalar, int Options, typename StorageIndex>
	void write_matlab(
		const std::string& file_name,
		Eigen::SparseMatrix<Scalar, Options, StorageIndex>& M)
	{
		std::ofstream writer(file_name);
		for (Eigen::Index outer = 0; outer < M.outerSize(); ++outer)
		{
			for (typename Eigen::SparseMatrix<Scalar, Options, StorageIndex>::InnerIterator iter(M, outer); iter; ++iter)
			{
				writer << (iter.row() + 1) << " " << (iter.col() + 1) << " ";
				if constexpr (na::is_complex_v<Scalar>)
				{
					writer << iter.value().real() << " " << iter.value().imag() << "\n";
				}
				else
				{
					writer << iter.value() << "\n";
				}
			}
		}
	}
}

#endif // !NA_IO_MATLAB_H
