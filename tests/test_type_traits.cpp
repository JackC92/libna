#include <complex>
#include <type_traits>
#include "gtest/gtest.h"
#include "na/type_traits/floating_point_scalar.h"
#include "na/type_traits/is_complex.h"
#include "na/type_traits/is_float_or_complex.h"

TEST(type_traits, is_complex)
{
	EXPECT_FALSE(na::is_complex_v<int>);
	EXPECT_FALSE(na::is_complex_v<float>);
	EXPECT_FALSE(na::is_complex_v<double>);
	EXPECT_FALSE(na::is_complex_v<long double>);
	EXPECT_TRUE(na::is_complex_v<std::complex<float>>);
	EXPECT_TRUE(na::is_complex_v<std::complex<double>>);
	EXPECT_TRUE(na::is_complex_v<std::complex<long double>>);
}

TEST(type_traits, is_floating_point_or_complex)
{
	EXPECT_FALSE(na::is_float_or_complex_v<int>);
	EXPECT_TRUE(na::is_float_or_complex_v<float>);
	EXPECT_TRUE(na::is_float_or_complex_v<double>);
	EXPECT_TRUE(na::is_float_or_complex_v<long double>);
	EXPECT_TRUE(na::is_float_or_complex_v<std::complex<float>>);
	EXPECT_TRUE(na::is_float_or_complex_v<std::complex<double>>);
	EXPECT_TRUE(na::is_float_or_complex_v<std::complex<long double>>);
}

TEST(type_traits, floating_point_scalar)
{
	bool float_test = std::is_same_v<na::floating_point_scalar<float>::value_type, float>;
	bool double_test = std::is_same_v<na::floating_point_scalar<double>::value_type, double>;
	bool long_double_test = std::is_same_v<na::floating_point_scalar<long double>::value_type, long double>;
	bool cmplx_float_test = std::is_same_v<na::floating_point_scalar<std::complex<float>>::value_type, float>;
	bool cmplx_double_test = std::is_same_v<na::floating_point_scalar<std::complex<double>>::value_type, double>;
	bool cmplx_long_double_test = std::is_same_v<na::floating_point_scalar<std::complex<long double>>::value_type, long double>;
	EXPECT_TRUE(float_test);
	EXPECT_TRUE(double_test);
	EXPECT_TRUE(long_double_test);
	EXPECT_TRUE(cmplx_float_test);
	EXPECT_TRUE(cmplx_double_test);
	EXPECT_TRUE(cmplx_long_double_test);
}
