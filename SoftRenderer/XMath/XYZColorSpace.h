// fitted CIE XYZ color matching functions by Nvidia,
// see "Simple Analytic Approximations to the CIE XYZ Color Matching Functions" @ Journal of Computer Graphics Techniques

#pragma once

#include <math.h>

namespace xmath
{
	namespace xyz1931
	{
		inline float X(float wavelen)
		{
			float dParam1 = (wavelen - 442.0f)*((wavelen < 442.0f) ? 0.0624f : 0.0374f);
			float dParam2 = (wavelen - 599.8f)*((wavelen < 599.8f) ? 0.0264f : 0.0323f);
			float dParam3 = (wavelen - 501.1f)*((wavelen < 501.1f) ? 0.0490f : 0.0382f);
			return 0.362f*expf(-0.5f*dParam1*dParam1) + 1.056f*expf(-0.5f*dParam2*dParam2) - 0.065f*expf(-0.5f*dParam3*dParam3);
		}

		inline float Y(float wavelen)
		{
			float dParam1 = (wavelen - 568.8f)*((wavelen < 568.8f) ? 0.0213f : 0.0247f);
			float dParam2 = (wavelen - 530.9f)*((wavelen < 530.9f) ? 0.0613f : 0.0322f);
			return 0.821f*expf(-0.5f*dParam1*dParam1) + 0.286f*expf(-0.5f*dParam2*dParam2);
		}

		inline float Z(float wavelen)
		{
			float dParam1 = (wavelen - 437.0f)*((wavelen < 437.0f) ? 0.0845f : 0.0278f);
			float dParam2 = (wavelen - 459.0f)*((wavelen < 459.0f) ? 0.0385f : 0.0725f);
			return 1.217f*expf(-0.5f*dParam1*dParam1) + 0.681f*expf(-0.5f*dParam2*dParam2);
		}
	}

	namespace xyz1964
	{
		inline float sqr(float f)
		{
			return f*f;
		}

		inline float X(float wavelen)
		{
			float val[4] = { 0.4f, 1014.0f, -0.02f, -570.f };
			float smallLobe = val[0] * expf(-1250 * sqr(logf((wavelen - val[3]) / val[1])));

			float val2[4] = { 1.13f, 234.f, -0.001345f, -1.799f };
			float bigLobe = val2[0] * expf(-val2[1] * sqr(logf((1338.0f - wavelen) / 743.5f)));

			return smallLobe + bigLobe;
		}

		inline float Y(float wavelen)
		{
			float val[4] = { 1.011f, 556.1f, 46.14f, 0.0f };
			return val[0] * expf(-0.5f* sqr((wavelen - val[1]) / val[2])) + val[3];
		}

		inline float Z(float wavelen)
		{
			float val[4] = { 2.06f, 180.4f, 0.125f, 266.0f };
			return val[0] * expf(-32.0f* sqr(logf((wavelen - val[3]) / val[1])));
		}
	}
}
