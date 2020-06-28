#pragma once
#include "Common.h"
#include "Vector.h"

#include <array>

namespace xmath
{
	template<typename _elem_ty, size_t _dimension>
	struct MatrixBase
	{
		using elem_ty = _elem_ty;
		using vec_ty = Vector<elem_ty, _dimension>;
		static const size_t dimension = _dimension;

		union
		{
			vec_ty row[dimension];
			elem_ty m[dimension][dimension];
		};

		MatrixBase(){}
	};

	template<>
	struct MatrixBase<float, 4>
	{
		using elem_ty = float;
		using vec_ty = Vector<elem_ty, 4>;
		static const size_t dimension = 4;

		union
		{
			vec_ty row[dimension];
			elem_ty m[dimension][dimension];

			struct
			{
				vec_ty xAxis;
				vec_ty yAxis;
				vec_ty zAxis;
				vec_ty position;
			};
		};

		MatrixBase() {}
	};
	


	// TODO support rvalueref construct and assignment
	template<typename _elem_ty = float, size_t _dimension = 4>
	struct Matrix :
		public MatrixBase<_elem_ty, _dimension>
	{
		using elem_ty = _elem_ty;
		using vec_ty = Vector<elem_ty, _dimension>;
		static const size_t dimension = _dimension;

		using super = MatrixBase<_elem_ty, _dimension>;

		/*
		using elem_ty = _elem_ty;
		using vec_ty = Vector<elem_ty, _dimension>;
		static const size_t dimension = _dimension;

		union
		{
			vec_ty row[dimension];
			elem_ty m[dimension][dimension];
		};
		*/

		Matrix()
		{
			for (size_t r = 0; r < dimension; ++r)
				for (size_t c = 0; c < dimension; ++c)
				{
					if (c == r)
						super::m[r][c] = elem_ty{1};
					else
						super::m[r][c] = elem_ty{0};
				}
		}

		Matrix(const Matrix& another)
		{
			for (size_t r = 0; r < dimension; ++r)
				super::row[r] = another.row[r];
		}

		Matrix(std::array<vec_ty, dimension> _rows)
		{
			for (size_t r = 0; r < dimension; ++r)
				super::row[r] = _rows[r];
		}

		vec_ty& operator [] (size_t i)
		{
			return super::row[i];
		}

		const vec_ty& operator [] (size_t i) const
		{
			return super::row[i];
		}

	#pragma region ArithmaticOperators

		Matrix & operator = (const Matrix & another)
		{
			for (size_t i = 0; i < dimension; ++i)
				super::row[i] = another.row[i];

			return *this;
		}

		Matrix operator * (const Matrix & another) const
		{
			Matrix ret;

			for (size_t r = 0; r < dimension; ++r)
			for (size_t c = 0; c < dimension; ++c)
			{
				ret[r][c] = 0;

				for (size_t i = 0; i < dimension; ++i)
				{
					ret[r][c] += super::row[r][i] * another[i][c];
				}
			}

			return ret;
		}

		Matrix operator * (elem_ty f) const
		{
			Matrix ret;

			for (size_t r = 0; r < dimension; ++r)
				for (size_t c = 0; c < dimension; ++c)
				{
					ret[r][c] = super::row[r][c] * f;
				}

			return ret;
		}

		Matrix operator / (elem_ty f) const
		{
			Matrix ret;

			for (size_t r = 0; r < dimension; ++r)
				for (size_t c = 0; c < dimension; ++c)
				{
					ret[r][c] = super::row[r][c] / f;
				}

			return ret;
		}

		Matrix operator + (const Matrix & another) const
		{
			Matrix ret;

			for (size_t r = 0; r < dimension; ++r)
				for (size_t c = 0; c < dimension; ++c)
				{
					ret[r][c] = super::row[r][c] + another.row[r][c];
				}

			return ret;
		}

		Matrix & operator *= (const Matrix & another)
		{
			*this = *this * another;

			return *this;
		}

		Matrix & operator += (const Matrix & another)
		{
			for (size_t r = 0; r < dimension; ++r)
				for (size_t c = 0; c < dimension; ++c)
				{
					super::row[r][c] += another.row[r][c];
				}

			return *this;
		}

	#pragma endregion

		void transpose()
		{
			for (size_t r = 0; r < dimension; ++r)
				for (size_t c = 0; c < r; ++c)
					swap(super::m[r][c], super::m[c][r]);
		}

	#pragma region MatrixCreation
		// TODO 这些矩阵生成函数放在这里合适吗？尤其是当 Matrix 类还是一个类模板的情况下

		static Matrix zeroMatrix()
		{
			Matrix ret;

			for (size_t i = 0; i < dimension; ++i)
			{

				ret[i][i] = 0;
			}


			return ret; // TODO，这个函数内部新建了一个矩阵，然后在返回的时候会进行一次拷贝，加上 move 以提高效率？不知道这种情况编译器是否会优化掉
		}

		static Matrix identityMatrix()
		{
			return Matrix();
		}

		static Matrix translationMatrix(float3 translation)
		{
			static_assert(dimension == 4, "xuxing::math::Matrix::translationMatrix: Only 4x4 matrices can have translation.");

			Matrix ret;

			ret.position.x = translation.x;
			ret.position.y = translation.y;
			ret.position.z = translation.z;
			ret.position.w = 1.f;

			return ret; // TODO，这个函数内部新建了一个矩阵，然后在返回的时候会进行一次拷贝，加上 move 以提高效率？不知道这种情况编译器是否会优化掉
		}

		static Matrix scalingMatrix(float3 scaling)
		{
			Matrix ret;

			for (size_t i = 0; i < std::min<size_t>(3, dimension); ++i)
			{
				ret[i][i] = scaling[i];
			}

			return ret; // TODO，这个函数内部新建了一个矩阵，然后在返回的时候会进行一次拷贝，加上 move 以提高效率？不知道这种情况编译器是否会优化掉
		}

		// create a 2D rotation matrix
		static Matrix rotationMatrix(float radian) 
		{
			static_assert(dimension == 2, "xuxing::math::Matrix::rotationMatrix(radian): Only 2x2 matrices can represent 2D rotation (which doesn't need to specify a rotation axis).");

			Matrix ret;

			ret[0][0] = cos(radian); ret[0][1] = sin(radian);
			ret[1][0] = -ret[0][1]; ret[1][1] = ret[0][0];

			return ret;
		}

		static Matrix rotationMatrix(float radian, float3 axis)
		{
			static_assert(dimension == 3 || dimension == 4, "xuxing::math::Matrix::rotationMatrix(radian, axis): Only 3x3 or 4x4 matrices can represent rotation around an axis.");

			Matrix ret;

			axis = normalize(axis);

			float x = axis.x;
			float y = axis.y;
			float z = axis.z;

			float xx = x * x;
			float yy = y * y;
			float zz = z * z;

			float xy = x * y;
			float xz = x * z;
			float yz = y * z;

			float sinr = sin(radian);
			float cosr = cos(radian);
			float one_minus_cosr = 1 - cosr;

			ret[0][0] = xx * one_minus_cosr + cosr;
			ret[0][1] = xy * one_minus_cosr + z * sinr;
			ret[0][2] = xz * one_minus_cosr - y * sinr;
			

			ret[1][0] = xy * one_minus_cosr - z * sinr;
			ret[1][1] = yy * one_minus_cosr + cosr;
			ret[1][2] = yz * one_minus_cosr + x * sinr;
			

			ret[2][0] = xz * one_minus_cosr + y * sinr;
			ret[2][1] = yz * one_minus_cosr - x * sinr;
			ret[2][2] = zz * one_minus_cosr + cosr;
			
			if (dimension == 4)
			{
				ret[0][3] = 0.f;
				ret[1][3] = 0.f;
				ret[2][3] = 0.f;
				ret[3][0] = 0.f;
				ret[3][1] = 0.f;
				ret[3][2] = 0.f;
				ret[3][3] = 1.f;
			}
			

			return ret;
		}

		static Matrix rotationMatrix(float radian, float3 axis, float3 center)
		{
			static_assert(dimension == 4, "xuxing::math::Matrix::rotationMatrix(radian, axis, center): Only 4x4 matrices can have translation.");

			Matrix ret;

			// 先平移到0
			ret = translationMatrix(-center);

			// 再旋转
			ret = mul(ret, rotationMatrix(radian, axis));

			// 再平移回 center 
			ret = mul(ret, translationMatrix(center));

			return ret;
		}

		// z is the camera's front direction
		static Matrix perspectiveProjectionMatrixLH(float fovy, float aspect, float zn, float zf)
		{
			static_assert(dimension == 4, "xuxing::math::Matrix::perspectiveProjectionMatrixLH: Only 4x4 matrices can be projection matrix.");

			auto h = 1.f / tan(fovy / 2.f);
			auto w = h / aspect;

			Matrix ret;

			ret[0][0] = w;
			ret[1][1] = h;
			ret[2][2] = zf / (zf - zn);
			ret[2][3] = 1.f;
			ret[3][2] = -zn * ret[2][2];// -zn * zf / (zf - zn);
			ret[3][3] = 0.f;

			return ret;
		}

		// -z is the camera's front direction
		static Matrix perspectiveProjectionMatrixRH(float fovy, float aspect, float zn, float zf)
		{
			static_assert(dimension == 4, "xuxing::math::Matrix::perspectiveProjectionMatrixRH: Only 4x4 matrices can be projection matrix.");

			auto h = 1.f / tan(fovy / 2.f);
			auto w = h / aspect;

			Matrix ret;

			ret[0][0] = w;
			ret[1][1] = h;
			ret[2][2] = zf / (zn - zf);
			ret[2][3] = -1.f;
			ret[3][2] = zn * ret[2][2];// zn * zf / (zn - zf);
			ret[3][3] = 0.f;

			return ret;
		}

		// z is the camera's front direction
		static Matrix orthoProjectionMatrixLH(float width, float height, float zn, float zf)
		{
			static_assert(dimension == 4, "xuxing::math::Matrix::orthoProjectionMatrixLH: Only 4x4 matrices can be projection matrix.");

			Matrix ret;

			ret[0][0] = 2.f / width;
			ret[1][1] = 2.f / height;
			ret[2][2] = 1.f / (zf - zn);
			ret[3][2] = -zn * ret[2][2];// -zn / (zf - zn);

			return ret;
		}

		// -z is the camera's front direction
		static Matrix orthoProjectionMatrixRH(float width, float height, float zn, float zf)
		{
			static_assert(dimension == 4, "xuxing::math::Matrix::orthoProjectionMatrixRH: Only 4x4 matrices can be projection matrix.");

			Matrix ret;

			ret[0][0] = 2.f / width;
			ret[1][1] = 2.f / height;
			ret[2][2] = 1.f / (zn - zf);
			ret[3][2] = zn * ret[2][2];// -zn / (zn - zf);

			return ret;
		}

	#pragma endregion
	};


	using float2x2 = Matrix<float, 2>;
	using float3x3 = Matrix<float, 3>;
	using float4x4 = Matrix<float, 4>;

	// vec is row vector
	template<typename vec_ty, typename mat_ty>
	vec_ty mul(const vec_ty & vec, const mat_ty & mat)
	{
		static_assert(vec_ty::dimension == mat_ty::dimension, "xuxing::math::mul: Dimension mismatch.");

		vec_ty ret;

		for (size_t c = 0; c < mat_ty::dimension; ++c)
		{
			for (size_t i = 0; i < mat_ty::dimension; ++i)
			{
				ret[c] += vec[i] * mat[i][c];
			}
		}

		return ret;
	}

	// TODO change name to "mul", use enable_if
	template<typename vec_ty, typename mat_ty>
	vec_ty mulColVector(const mat_ty & mat, const vec_ty & vec)
	{
		static_assert(vec_ty::dimension == mat_ty::dimension, "xuxing::math::mul: Dimension mismatch.");

		vec_ty ret;

		for (size_t r = 0; r < mat_ty::dimension; ++r)
		{
			for (size_t i = 0; i < mat_ty::dimension; ++i)
			{
				ret[r] += vec[i] * mat[r][i];
			}
		}

		return ret;
	}

	template<typename mat_ty>
	mat_ty mul(const mat_ty & l, const mat_ty & r)
	{
		return l * r;
	}

#pragma region ExperimentalBinaryArithmaticOperators

	inline float4 operator * (const float4& v, const float4x4& mat)
	{
		return mul(v, mat);
	}

#pragma endregion

	template<typename mat_ty>
	mat_ty transpose(mat_ty mat)
	{
		for (size_t r = 0; r < mat_ty::dimension; ++r)
		for (size_t c = 0; c < r; ++c)
			swap(mat[r][c], mat[c][r]);

		return mat;
	}

	// 实际上经过我的证明，即使是错切变换的逆也是它的转制？
	// 但是这样就和何文峰说的“刚性变换”以及和维基说的“正交变换”冲突了，这个应该继续研究一下
	// can only inverse a rigid transform???
	template<typename mat_ty>
	mat_ty fastInverse(const mat_ty & mat)
	{
		static_assert(mat_ty::dimension == 4, "xuxing::math::fastInverse: Only 4x4 matrices is supported by this method.");

		mat_ty ret({
			mat_ty::vec_ty(mat[0][0], mat[1][0], mat[2][0], (typename mat_ty::elem_ty)0),
			mat_ty::vec_ty(mat[0][1], mat[1][1], mat[2][1], (typename mat_ty::elem_ty)0),
			mat_ty::vec_ty(mat[0][2], mat[1][2], mat[2][2], (typename mat_ty::elem_ty)0),
			mat_ty::vec_ty((typename mat_ty::elem_ty)0, (typename mat_ty::elem_ty)0, (typename mat_ty::elem_ty)0, (typename mat_ty::elem_ty)1)
		});

		auto translation = -mat[3];
		translation[3] = 0;

		translation = mul(translation, ret);

		ret[3] += translation;

		return ret;
	}

	template<typename elem_ty>
	float determinant(const Matrix<elem_ty, 2>& mat)
	{
		return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
	}

	template<typename mat_ty>
	float determinant(const mat_ty& mat)
	{
		// cofactor expansion

		float ret = 0;

		for (size_t col = 0; col < mat_ty::dimension; ++col)
		{
			ret += mat[0][col] * cofactor(mat, 0, col);
		}

		return ret;
	}

	template<typename mat_ty>
	float cofactor(const mat_ty& mat, size_t elem_row, size_t elem_col)
	{
		Matrix<typename mat_ty::elem_ty, mat_ty::dimension - 1> minorMatrix;

		// minor

		for (size_t col = 0; col < mat_ty::dimension; ++col)
			for (size_t row = 0; row < mat_ty::dimension; ++row)
			{
				if (row == elem_row || col == elem_col)
					continue;

				size_t _row = row;
				size_t _col = col;

				if (row > elem_row)
					_row--;

				if (col > elem_col)
					_col--;

				minorMatrix[_row][_col] = mat[row][col];
			}

		// cofactor

		return determinant(minorMatrix) * pow(-1.f, elem_row + elem_col);
	}

	template<typename mat_ty>
	mat_ty adjugate(const mat_ty & mat)
	{
		mat_ty adjugateMatrix;

		for (size_t row = 0; row < mat_ty::dimension; ++row)
			for (size_t col = 0; col < mat_ty::dimension; ++col)
			{
				adjugateMatrix[row][col] = cofactor(mat, col, row);
			}

		return adjugateMatrix;
	}

	template<typename mat_ty>
	mat_ty inverse(const mat_ty & mat, float* pDeterminant = nullptr)
	{
		auto det = determinant(mat);

		if (pDeterminant)
			*pDeterminant = det;

		if (det == 0)
			return mat_ty::identityMatrix();

		return adjugate(mat) / det;
	}

	// x-y-z intrinsic, just like 3ds max's default, euler angles in radian
	inline void rotationMatrixToEulerAngle(float4x4 matrix, float3& rotation)
	{
		// x-y-z intrinsic = z-y-x extrinsic
		// matrix = Z(c) * Y(b) * X(a) 
		// see https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles

		float a, b, c;

		// cos_b != 0
		if ((matrix.m[0][0] + matrix.m[1][0]) != 0)
		{
			// assume cos_b > 0, b in (-PI/2, PI/2)

			float sin_b = matrix.m[2][0];
			b = asin(sin_b);
			float cos_b = cos(b);

			a = atan2(-matrix.m[2][1], matrix.m[2][2]);
			c = atan2(-matrix.m[1][0], matrix.m[0][0]);

			float sin_a = sin(a);
			float cos_a = cos(a);

			float sin_c = sin(c);
			float cos_c = cos(c);

			// under the assumption that cos_b > 0, 
			// calc some matrix elements and compare them with the given matrix

			float _01 = cos_a * sin_c + sin_a * sin_b * cos_c;	float _02 = sin_a * sin_c - cos_a * sin_b * cos_c;
			float _11 = cos_a * cos_c - sin_a * sin_b * sin_c;	float _12 = sin_a * cos_c + cos_a * sin_b * sin_c;

			if (abs(_01 - matrix.m[0][1]) > 0.001 || abs(_02 - matrix.m[0][2]) > 0.001 ||
				abs(_11 - matrix.m[1][1]) > 0.001 || abs(_12 - matrix.m[1][2]) > 0.001)
			{
				// elements do not match, assumption is false, so:
				// cos_b < 0,  b in (-PI, -PI/2) U (PI/2, PI)

				float sin_b = -matrix.m[2][0];
				float b = asin(sin_b);

				if (sin_b > 0 && b < pi / 2) // b should in (PI/2, PI)
					b = pi - b;

				if (sin_b < 0 && b > -pi / 2) // b should in (-PI, -PI/2)
					b = -pi - b;

				float cos_b = cos(b);

				a = atan2(matrix.m[2][1], -matrix.m[2][2]);
				c = atan2(matrix.m[1][0], -matrix.m[0][0]);
			}
		}
		// cos_b = 0; 
		else
		{
			// we are in gimbal lock, 
			// only b has a unique solution

			float sin_b = -matrix.m[2][0];
			float b = asin(sin_b);

			if (sin_b > 0)
				b = pi / 2;
			else
				b = -pi / 2;

			// a and c have no unique solution
			// assume c == 0
			c = 0;

			float sin_c = 0;
			float cos_c = 1;

			if (sin_b > 0) // sin_b == 1
			{
				a = atan2(matrix.m[0][1], matrix.m[1][1]);
			}
			else
			{
				a = atan2(matrix.m[0][2], matrix.m[1][2]);
			}
		}

		rotation.x = a;
		rotation.y = b;
		rotation.z = c;
	}

	// x-y-z intrinsic, just like 3ds max's default
	inline float4x4 eulerAngleToRotationMatrix(float3 rotation)
	{
		// intrinsic x-y-z = extrinsic z-y-x

		return 
			float4x4::rotationMatrix(rotation.z, float3(0, 0, 1)) * 
			float4x4::rotationMatrix(rotation.y, float3(0, 1, 0)) *
			float4x4::rotationMatrix(rotation.x, float3(1, 0, 0));
	}

	inline void decompose(float4x4 matrix, float3& position, float3& rotation, float3& scale)
	{
		// position

		position = matrix.position.toVector3();

		matrix.position = float4(0,0,0,1);

		// scale

		scale.x = length(matrix.xAxis);
		scale.y = length(matrix.yAxis);
		scale.z = length(matrix.zAxis);

		matrix.xAxis /= scale.x;
		matrix.yAxis /= scale.y;
		matrix.zAxis /= scale.z;

		// rotation

		rotationMatrixToEulerAngle(matrix, rotation);
	}

	inline float4x4 compose(float3 position, float3 rotation, float3 scale)
	{
		return float4x4::scalingMatrix(scale) * eulerAngleToRotationMatrix(rotation) * float4x4::translationMatrix(position);
	}


	// experimental helper functions

	inline float3 mulAsPosition(const float3& position, const float4x4& matrix)
	{
		return mul(float4(position, 1.f), matrix).toVector3();
	}

	inline float3 mulAsVector(const float3& vector, const float4x4& matrix)
	{
		return mul(float4(vector, 0.f), matrix).toVector3();
	}

	// experimental, column major matrix in memory
	template<typename _elem_ty = float, size_t _dimension = 4>
	struct ColMajorMatrix
	{
		using elem_ty = _elem_ty;
		using vec_ty = Vector<elem_ty, _dimension>;
		static const size_t dimension = _dimension;

		vec_ty col[dimension];

		ColMajorMatrix()
		{
			for (size_t r = 0; r < dimension; ++r)
				for (size_t c = 0; c < dimension; ++c)
				{
					if (c == r)
						col[r][c] = elem_ty{ 1 };
					else
						col[r][c] = elem_ty{ 0 };
				}
		}

		ColMajorMatrix(const Matrix<_elem_ty, _dimension>& rowMajor)
		{
			for (size_t r = 0; r < _dimension; ++r)
				for (size_t c = 0; c < _dimension; ++c)
				{
					col[c][r] = rowMajor.row[r][c];
				}
		}

		ColMajorMatrix& operator = (const Matrix<_elem_ty, _dimension>& rowMajor)
		{
			for (size_t r = 0; r < _dimension; ++r)
				for (size_t c = 0; c < _dimension; ++c)
				{
					col[c][r] = rowMajor.row[r][c];
				}

			return *this;
		}

		Matrix<_elem_ty, _dimension> rowMajor() const
		{
			Matrix<_elem_ty, _dimension> ret;

			for (size_t r = 0; r < _dimension; ++r)
				for (size_t c = 0; c < _dimension; ++c)
				{
					ret.row[r][c] = col[c][r];
				}

			return ret;
		}
	};

	using float2x2_col = ColMajorMatrix<float, 2>;
	using float3x3_col = ColMajorMatrix<float, 3>;
	using float4x4_col = ColMajorMatrix<float, 4>;
}
