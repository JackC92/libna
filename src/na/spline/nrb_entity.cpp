#include "na/spline/nrb_entity.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <numeric>
#include <vector>
#include "Eigen/Core"

namespace na
{
	namespace nurbs
	{
		NURBSEntity::NURBSEntity(
			const std::vector<Eigen::Index>& degree,
			const std::vector<Eigen::Index>& number,
			const std::vector<Eigen::VectorXd>& knots,
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const bool normalize)
		{
			assert(((number.size() == degree.size()) && (number.size() == knots.size())) && "NURBSEntity: degree, number and knots must have the same size");
			assert((coefs.rows() == std::accumulate(number.begin(), number.end(), Eigen::Index(1), std::multiplies<Eigen::Index>())) && "NURBSEntity: the number of rows of coefs must be compatible with number");
			assert(((coefs.cols() >= 1) && coefs.cols() <= 4) && "NURBSEntity: the number of columns of coefs must be between 1 and 4");

			Eigen::Index dim = coefs.cols();
			if (dim < 4)
			{
				m_coefs.setZero(std::accumulate(number.begin(), number.end(), Eigen::Index(1), std::multiplies<Eigen::Index>()), 4);
				m_coefs.leftCols(dim) = coefs;
				m_coefs.col(3).setConstant(1.0);
			}
			else
			{
				// coefs are given as homogeneous coordinates with the NURBS weights as the last column
				m_coefs = coefs;
			}

			m_number = number;
			for (std::size_t i = 0; i < m_number.size(); ++i)
			{
				m_degree.emplace_back(knots[i].size() - number[i] - 1);
				m_knots.emplace_back(knots[i]);
				std::sort<double*>(m_knots[i].data(), m_knots[i].data() + m_knots[i].size());
				if (normalize)
				{
					m_knots[i] = (m_knots[i].array() - m_knots[i](m_degree[i])) / (*(m_knots[i].end() - m_degree[i] - 1) - m_knots[i](m_degree[i]));
				}
				assert(((m_knots[i].array() == *m_knots[i].begin()).count() == m_degree[i] + 1) && "NURBSEntity: the multiplicity of the first entry must be degree + 1");
				assert(((m_knots[i].array() == *(m_knots[i].end() - 1)).count() == m_degree[i] + 1) && "NURBSEntity: the multiplicity of the last entry must be degree + 1");
			}
		}

		NURBSEntity::NURBSEntity(
			const Eigen::Index degree,
			const Eigen::Index number,
			const Eigen::Ref<const Eigen::VectorXd>& knots,
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const bool normalize)
			: NURBSEntity(std::vector<Eigen::Index>(1, degree), std::vector<Eigen::Index>(1, number), std::vector<Eigen::VectorXd>(1, knots), coefs)
		{

		}

		NURBSEntity::~NURBSEntity()
		{

		}

		std::vector<Eigen::Index>& NURBSEntity::degree()
		{
			return m_degree;
		}

		std::vector<Eigen::Index>& NURBSEntity::number()
		{
			return m_number;
		}

		std::vector<Eigen::VectorXd>& NURBSEntity::knots()
		{
			return m_knots;
		}

		Eigen::MatrixX4d& NURBSEntity::coefs()
		{
			return m_coefs;
		}

		const std::vector<Eigen::Index>& NURBSEntity::degree() const
		{
			return m_degree;
		}

		const std::vector<Eigen::Index>& NURBSEntity::number() const
		{
			return m_number;
		}

		const std::vector<Eigen::VectorXd>& NURBSEntity::knots() const
		{
			return m_knots;
		}

		const Eigen::MatrixX4d& NURBSEntity::coefs() const
		{
			return m_coefs;
		}
	}
}
