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
			const std::vector<Eigen::ArrayXd>& knots,
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const bool normalize)
		{
			Eigen::Index dim = coefs.cols();
			Eigen::Index ncp = std::accumulate(number.begin(), number.end(), static_cast<Eigen::Index>(1), std::multiplies<Eigen::Index>());

			assert(((degree.size() == number.size()) && (degree.size() == knots.size())) && "NURBSEntity: degree, number and knots must have the same size");
			assert((coefs.rows() == ncp) && "NURBSEntity: the number of rows of coefs must be compatible with number");
			assert(((dim >= 1) && dim <= 4) && "NURBSEntity: the number of columns of coefs must be between 1 and 4");

			if (dim < 4)
			{
				m_coefs.setZero(ncp, 4);
				m_coefs.leftCols(dim) = coefs;
				m_coefs.col(3).setConstant(1.0);
			}
			else
			{
				// coefs are given as homogeneous coordinates with the NURBS weights as the last column
				m_coefs = coefs;
			}

			m_degree = degree;
			m_number = number;
			m_knots = knots;
			for (std::size_t i = 0; i < m_number.size(); ++i)
			{
				assert((m_knots[i].size() == m_degree[i] + number[i] + 1) && "NURBSEntity: the number of control points must be compatible");
				std::sort<double*>(m_knots[i].data(), m_knots[i].data() + m_knots[i].size());
				if (normalize)
				{
					m_knots[i] = (m_knots[i] - m_knots[i].coeffRef(m_degree[i])) / (*(m_knots[i].end() - m_degree[i] - 1) - m_knots[i].coeffRef(m_degree[i]));
				}
			}
		}

		NURBSEntity::NURBSEntity(
			const Eigen::Index degree,
			const Eigen::Index number,
			const Eigen::Ref<const Eigen::ArrayXd>& knots,
			const Eigen::Ref<const Eigen::MatrixXd>& coefs,
			const bool normalize)
			: NURBSEntity(std::vector<Eigen::Index>({ degree }), std::vector<Eigen::Index>({ number }), std::vector<Eigen::ArrayXd>({ knots }), coefs)
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

		std::vector<Eigen::ArrayXd>& NURBSEntity::knots()
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

		const std::vector<Eigen::ArrayXd>& NURBSEntity::knots() const
		{
			return m_knots;
		}

		const Eigen::MatrixX4d& NURBSEntity::coefs() const
		{
			return m_coefs;
		}
	}
}
