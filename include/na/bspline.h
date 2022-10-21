#ifndef NA_BSPLINE_H
#define NA_BSPLINE_H
#include "Eigen/Core"

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
	}
}

#endif // !NA_BSPLINE_H
