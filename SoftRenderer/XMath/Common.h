#pragma once
#include <cmath>
#include <algorithm>

namespace xmath
{
	static const float pi = 3.141592653589793238462643383279502884197f;
	static const float e = 2.7182818284590452353602874713526624977572f;

	template<typename ty>
	inline ty radianToDegree(const ty & radian)
	{
		return radian / pi * 180.f;
	}

	template<typename ty>
	inline ty degreeToRadian(const ty & degree)
	{
		return degree / 180.f * pi;
	}

	template<typename ty>
	ty lerp(const ty & v0, const ty & v1, float t)
	{
		return v0 + t * (v1 - v0);
	}

	template<typename ty>
	ty clamp(const ty& v, const ty& _min, const ty& _max)
	{
		if (v < _min)
			return _min;
		else if (v > _max)
			return _max;
		else
			return v;
	}

	template<typename ty>
	ty saturate(const ty& v)
	{
		return clamp(v, ty(0), ty(1));
	}

	template<typename ty>
	ty sign(const ty& v)
	{
		if (v < 0)
			return ty(-1);
		else
			return ty(1);
	}

	template<typename ty>
	ty normalizePeriodical(ty v, ty period, ty minValue = 0)
	{
		v = v - minValue;

		period = abs(period);
		v = remainder(v, period);

		if (v < 0)
			v += period;
		
		return v + minValue;
	}
}
