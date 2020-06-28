#pragma once
#include <cmath>
#include <stdint.h>

namespace xmath
{
	template<typename _elem_ty, size_t _dimension>
	struct Vector
	{
		Vector()
		{
		//	static_assert(false, "xuxing::math::Vector: Type not supported, do your own specialization.");
		}
	};

	// specialization of 2 dimensional vector
	template<typename _elem_ty>
	struct Vector<_elem_ty, 2>
	{
		using elem_ty = _elem_ty;
		static const size_t dimension = 2;

		union
		{
			struct
			{
				elem_ty x, y;
			};
			struct
			{
				elem_ty r, g;
			};
			elem_ty v[2];
		};

		Vector() :
			x((elem_ty)0), y((elem_ty)0) {}

		Vector(const Vector & r) :
			x(r.x), y(r.y) {}

		Vector(elem_ty _x, elem_ty _y)
		{
			set(_x, _y);
		}

		void set(elem_ty _x, elem_ty _y)
		{
			x = (elem_ty)_x;
			y = (elem_ty)_y;
		}

		elem_ty * data() const
		{
			return &x;
		}

		size_t size() const // constexpr when c++11 avilable // number of elems 
		{
			return 2;
		}

		elem_ty & operator[](size_t i)
		{
			return v[i];
		}

		const elem_ty & operator[](size_t i) const
		{
			return v[i];
		}

	#pragma region Operators
		// assignment operators
		//

		Vector & operator = (const Vector r)
		{
			x = r.x;
			y = r.y;

			return *this;
		}

		Vector & operator += (const Vector & r)
		{
			x += r.x;
			y += r.y;

			return *this;
		}
		Vector & operator -= (const Vector & r)
		{
			x -= r.x;
			y -= r.y;

			return *this;
		}
		Vector & operator *= (const Vector & r)
		{
			x *= r.x;
			y *= r.y;

			return *this;
		}
		Vector & operator /= (const Vector & r)
		{
			x /= r.x;
			y /= r.y;

			return *this;
		}

		Vector & operator *= (elem_ty r)
		{
			x *= r;
			y *= r;

			return *this;
		}
		Vector & operator /= (elem_ty r)
		{
			x /= r;
			y /= r;

			return *this;
		}

		// unary operators
		//

		Vector operator + () const
		{
			return *this;
		}
		Vector operator - () const
		{
			return Vector(-x, -y);
		}

		// binary operators
		//

		Vector operator + (const Vector& r) const
		{
			return Vector(x + r.x, y + r.y);
		}
		Vector operator - (const Vector& r) const
		{
			return Vector(x - r.x, y - r.y);
		}
		Vector operator * (const Vector& r) const
		{
			return Vector(x * r.x, y * r.y);
		}
		Vector operator / (const Vector & r) const
		{
			return Vector(x / r.x, y / r.y);
		}

		Vector operator * (elem_ty r) const
		{
			return Vector(x * r, y * r);
		}
		Vector operator / (elem_ty r) const
		{
			return Vector(x / r, y / r);
		}

		bool operator == (const Vector& r) const
		{
			return (x == r.x) && (y == r.y);
		}
		bool operator != (const Vector& r) const
		{
			return x != r.x || y != r.y;
		}

	#pragma endregion
	};

	// specialization of 3 dimensional vector
	template<typename _elem_ty>
	struct Vector<_elem_ty, 3>
	{
		using elem_ty = _elem_ty;
		static const size_t dimension = 3;

		union
		{
			struct
			{
				elem_ty x, y, z;
			};
			struct
			{
				elem_ty r, g, b;
			};
			elem_ty v[3];
		};

		Vector() :
			x(elem_ty{ 0 }), y(elem_ty{ 0 }), z(elem_ty{ 0 }) {}

		Vector(const Vector & r) = default;
		Vector(Vector && r) = default;

		Vector(elem_ty all)
		{
			set(all, all, all);
		}

		Vector(elem_ty _x, elem_ty _y, elem_ty _z)
		{
			set(_x, _y, _z);
		}

		Vector(Vector<elem_ty, 2> vec2, elem_ty _z) :
			x(vec2.x), y(vec2.y), z(_z) {}

		void set(elem_ty _x, elem_ty _y, elem_ty _z)
		{
			x = _x;
			y = _y;
			z = _z;
		}

		Vector<elem_ty, 2> toVector2() const
		{
			return Vector<elem_ty, 2>(x, y);
		}

		elem_ty * data() const
		{
			return &x;
		}

		size_t size() const // constexpr when c++11 avilable // number of elems 
		{
			return 3;
		}

		elem_ty & operator[](size_t i)
		{
			return v[i];
		}

		const elem_ty & operator[](size_t i) const
		{
			return v[i];
		}

		// TODO change *= * to template function
	#pragma region Operators
		// assignment operators
		//

		Vector & operator = (const Vector& r) = default;


		Vector & operator += (const Vector & r)
		{
			x += r.x;
			y += r.y;
			z += r.z;

			return *this;
		}
		Vector & operator -= (const Vector & r)
		{
			x -= r.x;
			y -= r.y;
			z -= r.z;

			return *this;
		}
		Vector & operator *= (const Vector & r)
		{
			x *= r.x;
			y *= r.y;
			z *= r.z;

			return *this;
		}
		Vector & operator /= (const Vector & r)
		{
			x /= r.x;
			y /= r.y;
			z /= r.z;

			return *this;
		}

		Vector & operator *= (elem_ty r)
		{
			x *= r;
			y *= r;
			z *= r;

			return *this;
		}
		Vector & operator /= (elem_ty r)
		{
			x /= r;
			y /= r;
			z /= r;

			return *this;
		}

		// unary operators
		//

		Vector operator + () const
		{
			return *this;
		}
		Vector operator - () const
		{
			return Vector(-x, -y, -z);
		}

		// binary operators
		//

		Vector operator + (const Vector& r) const
		{
			return Vector(x + r.x, y + r.y, z + r.z);
		}
		Vector operator - (const Vector& r) const
		{
			return Vector(x - r.x, y - r.y, z - r.z);
		}
		Vector operator * (const Vector& r) const
		{
			return Vector(x * r.x, y * r.y, z * r.z);
		}
		Vector operator / (const Vector & r) const
		{
			return Vector(x / r.x, y / r.y, z / r.z);
		}

		Vector operator * (elem_ty r) const
		{
			return Vector(x * r, y * r, z * r);
		}
		Vector operator / (elem_ty r) const
		{
			return Vector(x / r, y / r, z / r);
		}

		bool operator == (const Vector& r) const
		{
			return (x == r.x) && (y == r.y) && (z == r.z);
		}
		bool operator != (const Vector& r) const
		{
			return x != r.x || y != r.y || z != r.z;
		}

	#pragma endregion
	};

	// specialization of 4 dimensional vector
	template<typename _elem_ty>
	struct Vector<_elem_ty, 4>
	{
		using elem_ty = _elem_ty;
		static const size_t dimension = 4;

		union
		{
			struct
			{
				elem_ty x, y, z, w;
			};
			struct
			{
				elem_ty r, g, b, a;
			};
			elem_ty v[4];
		};

		Vector() :
			x((elem_ty)0), y((elem_ty)0), z((elem_ty)0), w((elem_ty)0) {}

		Vector(const Vector & r) :
			x(r.x), y(r.y), z(r.z), w(r.w) {}

		Vector(elem_ty _x, elem_ty _y, elem_ty _z, elem_ty _w = 0)
		{
			set(_x, _y, _z, _w);
		}

		Vector(elem_ty all)
		{
			set(all, all, all, all);
		}

		template<typename r_elem_ty>
		Vector(Vector<r_elem_ty, 4> r):
			x(r.x), y(r.y), z(r.z), w(r.w) {}

		Vector(Vector<elem_ty, 3> vec3, elem_ty _w) :
			x(vec3.x), y(vec3.y), z(vec3.z), w(_w) {}

		template<typename ty>
		void set(ty _x, ty _y, ty _z, ty _w)
		{
			x = (elem_ty)_x;
			y = (elem_ty)_y;
			z = (elem_ty)_z;
			w = (elem_ty)_w;
		}

		Vector<elem_ty, 3> toVector3() const
		{
			return Vector<elem_ty, 3>(x, y, z);
		}

		// experimental
		operator Vector<elem_ty, 3>() const
		{
			return Vector<elem_ty, 3>(x, y, z);
		}


		elem_ty * data() const
		{
			return (elem_ty *)&x;
		}

		size_t size() const // constexpr when c++11 avilable // number of elems 
		{
			return 4;
		}

		elem_ty & operator[](size_t i)
		{
			return v[i];
		}

		const elem_ty & operator[](size_t i) const
		{
			return v[i];
		}

		// TODO change *= * to template function
		// TODO 这些应该大部分可以从特化里挪出去
	#pragma region Operators
		// assignment operators
		//

		Vector & operator = (const Vector r)
		{
			x = r.x;
			y = r.y;
			z = r.z;
			w = r.w;

			return *this;
		}

		Vector & operator += (const Vector & r)
		{
			x += r.x;
			y += r.y;
			z += r.z;
			w += r.w;

			return *this;
		}
		Vector & operator -= (const Vector & r)
		{
			x -= r.x;
			y -= r.y;
			z -= r.z;
			w -= r.w;

			return *this;
		}
		Vector & operator *= (const Vector & r)
		{
			x *= r.x;
			y *= r.y;
			z *= r.z;
			w *= r.w;

			return *this;
		}
		Vector & operator /= (const Vector & r)
		{
			x /= r.x;
			y /= r.y;
			z /= r.z;
			w /= r.w;

			return *this;
		}

		Vector & operator *= (elem_ty r)
		{
			x *= r;
			y *= r;
			z *= r;
			w *= r;

			return *this;
		}
		Vector & operator /= (elem_ty r)
		{
			x /= r;
			y /= r;
			z /= r;
			w /= r;

			return *this;
		}

		// unary operators
		//

		Vector operator + () const
		{
			return *this;
		}
		Vector operator - () const
		{
			return Vector(-x, -y, -z, -w);
		}

		// binary operators
		//

		Vector operator + (const Vector& r) const
		{
			return Vector(x + r.x, y + r.y, z + r.z, w + r.w);
		}
		Vector operator - (const Vector& r) const
		{
			return Vector(x - r.x, y - r.y, z - r.z, w - r.w);
		}
		Vector operator * (const Vector& r) const
		{
			return Vector(x * r.x, y * r.y, z * r.z, w * r.w);
		}
		Vector operator / (const Vector & r) const
		{
			return Vector(x / r.x, y / r.y, z / r.z, w / r.w);
		}

		Vector operator * (elem_ty r) const
		{
			return Vector(x * r, y * r, z * r, w * r);
		}
		Vector operator / (elem_ty r) const
		{
			return Vector(x / r, y / r, z / r, w / r);
		}

		bool operator == (const Vector& r) const
		{
			return (x == r.x) && (y == r.y) && (z == r.z) && (w == r.w);
		}
		bool operator != (const Vector& r) const
		{
			return x != r.x || y != r.y || z != r.z || w != r.w;
		}

	#pragma endregion
	};

#pragma region TypeDefs
	// using uint = unsigned int; // TODO this might conflict with Qt's uint, fuck those who don't use namespaces

	using float1 = float;
	using uint1 = uint32_t;
	using sint1 = int32_t;
	using ubyte1 = uint8_t;
	using sbyte1 = int8_t;

	using float2 = Vector<float, 2>;
	using uint2 = Vector<uint1, 2>;
	using sint2 = Vector<sint1, 2>;
	using ubyte2 = Vector<ubyte1, 2>;
	using sbyte2 = Vector<sbyte1, 2>;

	using float3 = Vector<float, 3>;
	using uint3 = Vector<uint1, 3>;
	using sint3 = Vector<sint1, 3>;
	using ubyte3 = Vector<ubyte1, 3>;
	using sbyte3 = Vector<sbyte1, 3>;

	using float4 = Vector<float, 4>;
	using uint4 = Vector<uint1, 4>;
	using sint4 = Vector<sint1, 4>;
	using ubyte4 = Vector<ubyte1, 4>;
	using sbyte4 = Vector<sbyte1, 4>;

#pragma endregion


	template<typename vec_ty, typename result_ty = float>
	result_ty dot(const vec_ty & l, const vec_ty & r)
	{
		result_ty result = 0;

		for (size_t i = 0; i < vec_ty::dimension; ++i)
			result += (result_ty)l.v[i] * (result_ty)r.v[i];

		return result;
	}

	// 3D vector cross, 4D vector is treated as 3D vector here
	template<typename vec_ty>
	vec_ty cross(const vec_ty & l, const vec_ty & r)
	{
		static_assert(vec_ty::dimension == 3 || vec_ty::dimension == 4, "xuxing::math::cross: Dimension not supported.");

		return vec_ty(l.y * r.z - l.z * r.y, l.z * r.x - l.x * r.z, l.x * r.y - l.y * r.x);
	}

	template<typename vec_ty, typename result_ty = float>
	result_ty length(const vec_ty & v)
	{
		result_ty len2 = 0;

		for (size_t i = 0; i < vec_ty::dimension; ++i)
			len2 += v[i] * v[i];

		return sqrt(len2);
	}

	template<typename vec_ty>
	vec_ty normalize(const vec_ty & v)
	{
		return v / length(v);
	}

	template<typename vec_ty>
	vec_ty reflect(const vec_ty & v, const vec_ty & n)
	{
		static_assert(vec_ty::dimension == 3 || vec_ty::dimension == 4, "xuxing::math::cross: Dimension not supported.");

		return v - 2 * dot(v, n) * n;
	}

	template<typename vec_ty>
	float distance(const vec_ty & l, const vec_ty & r)
	{
		return length(l-r);
	}


	template<typename vec_ty>
	vec_ty operator * (typename vec_ty::elem_ty l, const vec_ty & v)
	{
		return v * l;
	}
}
