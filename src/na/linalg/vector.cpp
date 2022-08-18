#include "na/linalg/vector.h"
#include <cmath>
#include "na/macros.h"

namespace na
{
	namespace linalg
	{
		EXPLICIT_SPEC float scaling_constants<float>::tau_min = std::pow(0.5f, 49.0f);
		EXPLICIT_SPEC float scaling_constants<float>::tau_max = std::pow(2.0f, 62.0f);
		EXPLICIT_SPEC float scaling_constants<float>::sig_min = std::pow(2.0f, 100.0f);
		EXPLICIT_SPEC float scaling_constants<float>::sig_max = std::pow(0.5f, 66.0f);

		EXPLICIT_SPEC double scaling_constants<double>::tau_min = std::pow(0.5, 482.0);
		EXPLICIT_SPEC double scaling_constants<double>::tau_max = std::pow(2.0, 510.0);
		EXPLICIT_SPEC double scaling_constants<double>::sig_min = std::pow(2.0, 592.0);
		EXPLICIT_SPEC double scaling_constants<double>::sig_max = std::pow(0.5, 514.0);
	}
}
