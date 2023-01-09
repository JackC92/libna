#if !defined(NA_MATLAB_CORE_H) && defined(NA_USE_MATLAB)
#define NA_MATLAB_CORE_H
#include <memory>
#include <string>
#include "MatlabDataArray.hpp"
#include "MatlabEngine.hpp"

namespace na
{
	namespace matlab
	{
		class MATLABSession
		{
		public:
			MATLABSession();
			MATLABSession(const std::u16string& name);
			~MATLABSession();

			std::unique_ptr<::matlab::engine::MATLABEngine>& get_engine();
			::matlab::data::Reference<::matlab::data::Struct> get_workspace();
			void set_workspace();
			void update_workspace();
			
			void update_variable(const std::string& name, const ::matlab::data::Array& value);

		private:
			void initialize();

			std::unique_ptr<::matlab::engine::MATLABEngine> m_engine;
			::matlab::data::Array m_workspace;
		};
	}
}

#endif // !defined(NA_MATLAB_CORE_H) && defined(NA_USE_MATLAB)
