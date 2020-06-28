#pragma once

#include "resource.h"
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <atlimage.h>

class Vector3;
class Vector4;
class Matrix44;


COLORREF LerpColor(COLORREF color0, COLORREF color1, float k);
float LerpFloat(float f0, float f1, float k);
Vector3 ColorToVector3(COLORREF c);

class Vector3
{
public:
	float x;
	float y;
	float z;

	Vector3(float x, float y, float z);
	Vector3(Vector4 v);

	static Vector3 Lerp(Vector3 v0, Vector3 v1, float k);
	static Vector3 HadamardProduct(Vector3 v0, Vector3 v1);
	static Vector3 Reflect(Vector3 v, Vector3 n);

	Vector3 operator+(const Vector3 &v);
	Vector3 operator-(const Vector3 &v);
	Vector3 operator-();
	float operator*(const Vector3 &v);
	Vector3 operator*(const float &f);
	Vector3 operator^(const Vector3 &v);
	float GetLength();
	Vector3 Normalize();

};

class Vector4
{
public:
	float x;
	float y;
	float z;
	float w;

	Vector4(float x, float y, float z, float w);
	Vector4(Vector3 v, float w);

	Vector4 operator+(const Vector4 &v);
	Vector4 operator-(const Vector4 &v);
	float operator*(const Vector4 &v);
	Vector4 operator*(const Matrix44 &v);
	Vector4 operator*(const float &f);
	Vector4 operator^(const Vector4 &v);
	float GetLength();
	Vector4 Normalize();

};

class Matrix44
{
public:
	float m[4][4];

	Matrix44();

	static Matrix44 Identity();

	Matrix44 operator*(const Matrix44 &v);
	Vector4 operator*(const Vector4 &v);
	Matrix44 Transpose();
	Matrix44 InverseDiagonal();

};


class Vector2
{
public:
	float x;
	float y;
	Vector2(float x, float y);

	static Vector2 Lerp(Vector2 v0, Vector2 v1, float k);

	Vector2 operator+(const Vector2 &v);
	Vector2 operator-(const Vector2 &v);
	float operator*(const Vector2 &v);
	Vector2 operator*(const float &f);
};

typedef Vector2 TexCoord;

class Vertex
{
public:
	Vector3 pos;
	Vector3 wPos;
	Vector3 normal;
	TexCoord texcoord;
	COLORREF color;
	float rhw;

	Vertex(Vector3 pos, TexCoord texcoord, Vector3 normal, COLORREF color);

	static Vertex Lerp(Vertex v0, Vertex v1, float k);
	//Vertex(Vector3 pos, TexCoord texcoord, Vector3 normal);
	//Vertex();
};

class Triangle
{
public:
	Vertex v[3];

	Triangle(Vertex v0, Vertex v1, Vertex v2);

};

class Camera {
public:
	Vector3 pos;
	Vector3 lookAt;
	Vector3 up;

	Camera(Vector3 pos, Vector3 lookAt, Vector3 up);
};

class Mesh
{
public:
	std::vector<Vector3> v;
	std::vector<TexCoord> vt;
	std::vector<Vector3> vn;
	std::vector<Triangle> f;

	void LoadOBJ(const char* filename);
	void SetValue(const std::string tag, std::vector<std::string> temp);
	void ParseFace(std::string s0, std::string s1, std::string s2);
	void SplitString(std::vector<std::string>& v, const std::string& s, const std::string& d);
};

class Transform
{
public:
	Matrix44 translate;
	Matrix44 rotation;
	Matrix44 scale;

	Camera camera;

	Matrix44 modelMatrix;
	Matrix44 viewMatrix;
	Matrix44 projectionMatrix;
	Matrix44 MVPMatrix;
	Matrix44 normalMatrix;

	Matrix44 viewportMatrix;

	Transform();
	void SetTranslate(float x, float y, float z);
	void SetRotation(float x, float y, float z);
	void SetScale(float x, float y, float z);
	void UpdateModelMatrix();
	void SetCamera(Vector3 pos, Vector3 lookAt, Vector3 up);
	void UpdateViewMatrix();
	void SetprojectionMatrix(float fovy, float w, float h, float zn, float zf);
	void UpdateMVP();
	void UpdateNormalMatrix();
	void SetViewportMatrix(float x, float y, float w, float h, float zmin, float zmax);
	void Update();
	Vector4 PerspectiveDivision(Vector4 v);
	Vertex Apply(Vertex v);
};





class RenderBuffer
{
public:
	int width;
	int height;
	std::vector<float> zBuffer;
	std::vector<COLORREF> frameBuffer;

	RenderBuffer(int w = 800, int h = 600);

	static TexCoord LerpTexCoord(Vertex v0, Vertex v1, float k);

	void Resize(int w, int h);
	void LineDDA(Vertex v0, Vertex v1);
	void LineBresenham(Vertex v0, Vertex v1);
	void FillTriangle(Vertex v0, Vertex v1, Vertex v2);
	void ScanLine(int y, int x0, int x1, Vertex v0, Vertex v1);
	void ScanLineTex(int y, int x0, int x1, Vertex v0, Vertex v1);
	void FillPixel(int x, int y, float z, COLORREF c);
	void DrawFrame(HDC hdc);
};

class DirectionalLight
{
public:
	Vector3 direction;
	Vector3 color;

	DirectionalLight();
	Vector3 Apply(Vertex v, Camera camera);

};


