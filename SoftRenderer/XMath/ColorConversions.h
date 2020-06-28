#pragma once

#include "XYZColorSpace.h"
#include "Vector.h"
#include "Matrix.h"

namespace xmath
{
	namespace xyz
	{
		inline float3 to_sRGB_linear(float3 xyz)
		{
			float3x3 M
			{ {
				float3{3.2404542f, -1.5371385f, -0.4985314f},
				float3{-0.9692660f, 1.8760108f, 0.0415560f},
				float3{0.0556434f, -0.2040259f, 1.0572252f}
			} };

			float3 sRGB_linear = mulColVector(M, xyz);
			sRGB_linear.x = std::max(sRGB_linear.x, 0.f);
			sRGB_linear.y = std::max(sRGB_linear.y, 0.f);
			sRGB_linear.z = std::max(sRGB_linear.z, 0.f);

			return sRGB_linear;
		}

		inline float3 from_sRGB_linear(float3 rgb)
		{
			float3x3 M
			{ {
				float3{0.4124564f, 0.3575761f, 0.1804375f},
				float3{0.2126729f, 0.7151522f, 0.0721750f},
				float3{0.0193339f, 0.1191920f, 0.9503041f}
			} };

			return mulColVector(M, rgb);
		}
	}
}
