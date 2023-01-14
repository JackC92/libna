#include "na/spline/bzr_evaluate.h"
#include <cmath>
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "na/special/gamma.h"

namespace na
{
	namespace bezier
	{
		Eigen::Vector3d bzr_evaluate(
			const Eigen::Ref<const Eigen::MatrixX3d>& coefs,
			const double u)
		{
			Eigen::Vector3d vec = Eigen::Vector3d::Zero();
			for (Eigen::Index i = 0; i < coefs.rows(); ++i)
			{
				vec += (na::gamma::binomial_coefficient(coefs.rows() - 1, i) * std::pow(u, i) * std::pow(1.0 - u, coefs.rows() - i - 1)) * coefs.row(i);
			}
			return vec;
		}
		
		Eigen::Vector3d bzr_evaluate(
			const Eigen::Ref<const Eigen::MatrixX4d>& coefs,
			const double u)
		{
			Eigen::Vector4d vec = Eigen::Vector4d::Zero();
			for (Eigen::Index i = 0; i < coefs.rows(); ++i)
			{
				vec += (na::gamma::binomial_coefficient(coefs.rows() - 1, i) * std::pow(u, i) * std::pow(1.0 - u, coefs.rows() - i - 1)) * coefs.row(i);
			}
			return vec.hnormalized();
		}
	}
}
