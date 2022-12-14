#include "na/special/bessel.h"
#include <cassert>
#include <cmath>

namespace na
{
	namespace bessel
	{
		double besselj0_zero(const int k)
		{
			// Computed using N[BesselJZero[0, k], 20] in Mathematica
			static double j0zero[20] =
			{
				2.4048255576957727686,
				5.5200781102863106496,
				8.6537279129110122170,
				11.791534439014281614,
				14.930917708487785948,
				18.071063967910922543,
				21.211636629879258959,
				24.352471530749302737,
				27.493479132040254796,
				30.634606468431975118,
				33.775820213573568684,
				36.917098353664043980,
				40.058425764628239295,
				43.199791713176730358,
				46.341188371661814019,
				49.482609897397817174,
				52.624051841114996029,
				55.765510755019979312,
				58.906983926080942133,
				62.048469190227169883
			};

			assert((k >= 1) && "besselj0_zero: k must be a positive integer");

			if (k > 20)
			{
				// For k > 20, the zeros are computed using McMahon's asymptotic expansions for large zeros (equation (10.21.19) in the NIST Handbook of Mathematical Functions).
				const double ak = M_PI * (k - 0.25);
				const double akinv = 1.0 / ak;
				const double akinv2 = akinv * akinv;
				return ak + akinv * (0.125 + akinv2 * (-8.0729166666666666667e-2 + akinv2 * (2.4602864583333333333e-1 + akinv2 * (-1.8244387672061010974 + akinv2 * (25.336414797343905010 + akinv2 * (-567.64441213518338114 + akinv2 * (18690.476528232065383 + akinv2 * (-8.4935358029914876992e5 + 5.0922546240222676950e7 * akinv2))))))));
			}
			else
			{
				return j0zero[k - 1];
			}
		}

		double besselj1_squared(const int k)
		{
			// Computed using N[BesselJ[1, BesselJZero[0, k]], 20] in Mathematica
			static double j1squared[20] =
			{
				2.6951412394191692614e-1,
				1.1578013858220369581e-1,
				7.3686351136408215141e-2,
				5.4037573198116282042e-2,
				4.2661429017243091266e-2,
				3.5242103490996101359e-2,
				3.0021070103054672675e-2,
				2.6147391495308088590e-2,
				2.3159121824691392265e-2,
				2.0783829122267857604e-2,
				1.8850450669317667816e-2,
				1.7246157569665008300e-2,
				1.5893518105923597803e-2,
				1.4737626096472189590e-2,
				1.3738465145387117918e-2,
				1.2866181737615132879e-2,
				1.2098051548626797547e-2,
				1.1416471224491608517e-2,
				1.0807592791180204012e-2,
				1.0260372926280762811e-2
			};

			assert((k >= 1) && "besselj1_squared: k must be a positive integer");

			if (k > 20)
			{
				const double x = 1.0 / (k - 0.25);
				const double x2 = x * x;
				return x * (2.0264236728467555082e-1 + x2 * x2 * (-3.0338042971129032619e-4 + x2 * (1.9892436424596935947e-4 + x2 * (-2.2896990277211174141e-4 + x2 * (4.3371071913074645378e-4 + x2 * (-1.2363234972717541473e-3 + x2 * (4.9610142326888310287e-3 + x2 * (-2.6683739370232375770e-2 + 1.8539539820634562871e-1 * x2))))))));
			}
			else
			{
				return j1squared[k - 1];
			}
		}
	}
}