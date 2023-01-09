#if !defined(NA_MATLAB_COMMANDS_H) && defined(NA_USE_MATLAB)
#define NA_MATLAB_COMMANDS_H
#include <cassert>
#include <string>
#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "na/core/string.h"
#include "na/macros.h"
#include "na/matlab/adaptor.h"
#include "na/matlab/core.h"

namespace na
{
	namespace matlab
	{
		template <typename DerivedIn, typename DerivedV, typename DerivedD>
		void eig(
			MATLABSession& session,
			const Eigen::DenseBase<DerivedIn>& mat,
			Eigen::PlainObjectBase<DerivedV>& V,
			Eigen::DiagonalBase<DerivedD>& D)
		{
			session.update_variable("temp", to_matlab(mat));

			std::string cmd = na::format("[%s.V, %s.D] = eig(%s.temp)", MATLAB_WORKSPACE_NAME, MATLAB_WORKSPACE_NAME, MATLAB_WORKSPACE_NAME);
			session.get_engine()->eval(std::u16string(cmd.cbegin(), cmd.cend()));
			session.update_workspace();

			assert((session.get_workspace()["V"].getType() == ::matlab::data::GetArrayType<typename DerivedD::Scalar>::type) && "eig: V is not of the correct numeric type");
			assert((session.get_workspace()["D"].getType() == ::matlab::data::GetArrayType<typename DerivedD::Scalar>::type) && "eig: D is not of the correct numeric type");

			::matlab::data::TypedArrayRef<typename DerivedV::Scalar> Vref = session.get_workspace()["V"];
			::matlab::data::TypedArrayRef<typename DerivedD::Scalar> Dref = session.get_workspace()["D"];
			to_eigen<DerivedV>(Vref, V);
			to_eigen<DerivedD>(Dref, D);
		}

		template <typename Derived>
		void eig(
			MATLABSession& session,
			const Eigen::DenseBase<Derived>& mat,
			const char* V,
			const char* D)
		{
			session.update_variable("temp", to_matlab(mat));

			std::string cmd = na::format("[%s, %s] = eig(%s.temp)", V, D, MATLAB_WORKSPACE_NAME);
			session.get_engine()->eval(std::u16string(cmd.cbegin(), cmd.cend()));
			session.update_workspace();
		}

		template <typename Scalar, int Options, typename StorageIndex>
		void spy(
			MATLABSession& session,
			const Eigen::SparseMatrix<Scalar, Options, StorageIndex>& mat)
		{
			session.update_variable("temp", to_matlab(mat));

			std::string cmd = na::format("spy(%s.temp)", MATLAB_WORKSPACE_NAME);
			session.get_engine()->eval(std::u16string(cmd.cbegin(), cmd.cend()));
		}

		template <typename Scalar, int Options, typename StorageIndex>
		void spy(
			MATLABSession& session,
			const std::u16string& name,
			const Eigen::SparseMatrix<Scalar, Options, StorageIndex>& mat)
		{
			session.get_engine()->setVariable(name, to_matlab(mat));
			session.get_engine()->eval(u"spy(" + name + u")");
		}
	}
}

#endif // !defined(NA_MATLAB_COMMANDS_H) && defined(NA_USE_MATLAB)
