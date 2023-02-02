#include "na/matlab/core.h"
#ifdef NA_MATLAB_CORE_H
#include <memory>
#include <string>
#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"
#include "na/macros.h"

namespace na
{
	namespace matlab
	{
		MATLABSession::MATLABSession()
		{
			m_engine = ::matlab::engine::connectMATLAB();
			initialize();
		}

		MATLABSession::MATLABSession(const std::u16string& name)
		{
			m_engine = ::matlab::engine::connectMATLAB(name);
			initialize();
		}

		MATLABSession::~MATLABSession()
		{
			//m_engine->eval(std::u16string(u"clear ") + MATLAB_WORKSPACE_NAME);
		}

		std::unique_ptr<::matlab::engine::MATLABEngine>& MATLABSession::get_engine()
		{
			return m_engine;
		}

		::matlab::data::Reference<::matlab::data::Struct> MATLABSession::get_workspace()
		{
			return (*static_cast<::matlab::data::StructArray*>(&m_workspace))[0].operator ::matlab::data::Reference<::matlab::data::Struct>();
		}

		void MATLABSession::set_workspace()
		{
			m_engine->setVariable(MATLAB_WORKSPACE_NAME, m_workspace);
		}

		void MATLABSession::update_workspace()
		{
			m_workspace = m_engine->getVariable(MATLAB_WORKSPACE_NAME);
		}

		void MATLABSession::update_variable(const std::string& name, const ::matlab::data::Array& value)
		{
			// This only updates the Struct on the C++ side, not the actual variable in MATLAB
			get_workspace()[name] = value;
			// This will update the actual variable in MATLAB.
			set_workspace();
		}

		void MATLABSession::initialize()
		{
			::matlab::data::ArrayFactory factory;
			::matlab::data::StructArray arr = factory.createStructArray({ 1, 1 }, { "temp" });
			arr[0]["temp"] = factory.createArray<double>({ 1, 1 }, { 0.0 });
			m_engine->setVariable(MATLAB_WORKSPACE_NAME, arr);
			m_workspace = m_engine->getVariable(MATLAB_WORKSPACE_NAME);
		}
	}
}

#endif // NA_MATLAB_CORE_H
