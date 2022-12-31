#ifndef NA_SPLINE_NRB_ENTITY_H
#define NA_SPLINE_NRB_ENTITY_H
#include <vector>
#include "Eigen/Core"

namespace na
{
	namespace nurbs
	{
		class NURBSEntity
		{
		public:
			NURBSEntity(
				const std::vector<Eigen::Index>& degree,
				const std::vector<Eigen::Index>& number,
				const std::vector<Eigen::VectorXd>& knots,
				Eigen::Ref<const Eigen::MatrixXd> coefs,
				const bool normalize = false);
			NURBSEntity(
				const Eigen::Index degree,
				const Eigen::Index number,
				Eigen::Ref<const Eigen::VectorXd> knots,
				Eigen::Ref<const Eigen::MatrixXd> coefs,
				const bool normalize = false);
			~NURBSEntity();

			std::vector<Eigen::Index>& degree();
			std::vector<Eigen::Index>& number();
			std::vector<Eigen::VectorXd>& knots();
			Eigen::MatrixX4d& coefs();
			const std::vector<Eigen::Index>& degree() const;
			const std::vector<Eigen::Index>& number() const;
			const std::vector<Eigen::VectorXd>& knots() const;
			const Eigen::MatrixX4d& coefs() const;

		private:
			std::vector<Eigen::Index> m_degree;
			std::vector<Eigen::Index> m_number;
			std::vector<Eigen::VectorXd> m_knots;
			Eigen::MatrixX4d m_coefs;
		};
	}
}

#endif // !NA_SPLINE_NRB_ENTITY_H
