#ifndef NA_TYPE_TRAITS_MATLAB_TRAITS_H
#define NA_TYPE_TRAITS_MATLAB_TRAITS_H
#include <complex>
#include "na/macros.h"

namespace na
{
	namespace matlab
	{
		template <typename T>
		struct sparse_buffer_traits {};

		EXPSPEC struct sparse_buffer_traits<bool>
		{
			using data_type = bool;
		};

		EXPSPEC struct sparse_buffer_traits<int>
		{
			using data_type = double;
		};

		EXPSPEC struct sparse_buffer_traits<float>
		{
			using data_type = double;
		};

		EXPSPEC struct sparse_buffer_traits<double>
		{
			using data_type = double;
		};

		EXPSPEC struct sparse_buffer_traits<long double>
		{
			using data_type = double;
		};

		template <typename T>
		struct sparse_buffer_traits<std::complex<T>>
		{
			using data_type = std::complex<double>;
		};
	}
}

#endif // !NA_TYPE_TRAITS_MATLAB_TRAITS_H
