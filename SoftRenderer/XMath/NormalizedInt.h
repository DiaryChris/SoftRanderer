#pragma once
#include <cmath>
#include <limits>
#include <type_traits>
#include "Vector.h"

#undef min
#undef max

namespace xmath
{
	using namespace std;

	// for signed int_ty, [int_ty min, int_ty max] => [-1.f, 1.f]
	// for unsigned int_ty, [0, int_ty max] => [0.f, 1.f]
	template<typename int_ty>
	struct NormalizedInt
	{
		NormalizedInt(){}

		NormalizedInt(const NormalizedInt<int_ty>& r):
			value(r.value) {}

		template<typename ty, enable_if_t<is_floating_point<ty>::value, int> = 0>
		NormalizedInt(ty n)
		{
			fromFloat(n);
		}
		
		template<typename ty, enable_if_t<is_integral<ty>::value, int> = 0>
		NormalizedInt(ty n)
		{
			value = n;
		}

		void fromFloat(float f)
		{
			if (is_signed<int_ty>::value)
				f = clamp(f, -1.f, 1.f);
			else
				f = clamp(f, 0.f, 1.f);

			value = numeric_limits<int_ty>::max() * f;
		}

		float toFloat() const
		{
			return (float)value / (float)numeric_limits<int_ty>::max();
		}

		int_ty value = 0;
	};

	template<typename int_ty>
	using norm = NormalizedInt<int_ty>; // same name as ms direct x

	template<typename int_ty, size_t dim>
	using norm_vec = Vector<norm<int_ty>, dim>;

	template<typename int_ty>
	using norm3 = Vector<norm<int_ty>, 3>;

	template<typename int_ty>
	using norm4 = Vector<norm<int_ty>, 4>;
}
