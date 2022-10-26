#ifndef NA_BSPLINE_H
#define NA_BSPLINE_H
#include <cmath>
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Geometry"

namespace na
{
	namespace bspline
	{
		Eigen::Vector3d evaluate(
			const Eigen::MatrixXd& q,
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s);

		void evaluate_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s,
			int& mu,
			Eigen::VectorXd& b);

		Eigen::VectorXd evaluate_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s);

		Eigen::MatrixXd evaluate_basis(
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const Eigen::VectorXd& s);

		template <typename Node, typename MemPtr>
		Eigen::Vector3d evaluate(
			const std::vector<Node>& nodes,
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s,
			MemPtr func)
		{
			static_assert(std::is_member_function_pointer_v<MemPtr>, "func must be a non-static member function");

			auto fn = std::mem_fn(func);

			Eigen::VectorXd b;
			int mu;
			evaluate_basis(T, deg, m, s, mu, b);

			Eigen::Vector3d val(0.0, 0.0, 0.0);
			for (int i = 0; i < deg + 1; ++i)
			{
				val += fn(nodes[mu - deg + i]) * b(i);
			}
			return val;
		}

		template <typename Node, typename MemPtr>
		Eigen::Vector3d evaluate(
			const std::vector<std::shared_ptr<Node>>& nodes,
			const Eigen::VectorXd& T,
			const int deg,
			const int m,
			const double s,
			MemPtr func)
		{
			static_assert(std::is_member_function_pointer_v<MemPtr>, "func must be a non-static member function");

			auto fn = std::mem_fn(func);

			Eigen::VectorXd b;
			int mu;
			evaluate_basis(T, deg, m, s, mu, b);

			Eigen::Vector3d val(0.0, 0.0, 0.0);
			for (int i = 0; i < deg + 1; ++i)
			{
				val += fn(nodes[mu - deg + i].get()) * b(i);
			}
			return val;
		}

		Eigen::Vector3d curvature_binormal(
			const Eigen::MatrixXd& q,
			const Eigen::VectorXd& T,
			const int deg,
			const double s);

		template <typename Node, typename MemPtr>
		Eigen::Vector3d curvature_binormal(
			const std::vector<Node>& nodes,
			const Eigen::VectorXd& T,
			const int deg,
			const double s,
			MemPtr func)
		{
			Eigen::Vector3d v1 = evaluate(nodes, T, deg, 1, s, func);
			Eigen::Vector3d v2 = evaluate(nodes, T, deg, 2, s, func);
			return v1.cross(v2) * std::pow(v1.squaredNorm(), -1.5);
		}
	}
}

#endif // !NA_BSPLINE_H
