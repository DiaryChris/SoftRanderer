// SoftRenderer.cpp : 定义应用程序的入口点。
//



#include "stdafx.h"
#include "math.h"
#include "SoftRenderer.h"

#define MAX_LOADSTRING 100

// 全局变量:
HINSTANCE hInst;                                // 当前实例
WCHAR szTitle[MAX_LOADSTRING];                  // 标题栏文本
WCHAR szWindowClass[MAX_LOADSTRING];            // 主窗口类名


const char* filename = "sphere.obj";
const wchar_t* texFilename = _T("1bw.bmp");
bool isTex = true;
bool isPerspectiveCorrectOn = true;
CImage tex;
int texWidth;
int texHeight;
Transform trans;
Mesh mesh;
RenderBuffer rb;
DirectionalLight light;
Vector3 ambientColor = Vector3(1, 1, 1);
Vector3 lightDirection = Vector3(1, -1, -1);
Vector3 lightColor = Vector3(0.8f, 1, 1);
float ambientStrength = 0.1f;
float shininess = 16;
float rotationX = 130;
float rotationY = 200;
float rotationStep = 2;
float cameraZ = 3;
float cameraX = 0;
float cameraStep = 0.1f;
float cameraLookAtX = 0;
float cameraLookAtStep = 0.1f;


// 此代码模块中包含的函数的前向声明:
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);




#pragma region Tools
//颜色线性插值
COLORREF LerpColor(COLORREF color0, COLORREF color1, float k)
{

	float colorR = GetRValue(color1) * k + GetRValue(color0) * (1 - k);
	float colorG = GetGValue(color1) * k + GetGValue(color0) * (1 - k);
	float colorB = GetBValue(color1) * k + GetBValue(color0) * (1 - k);
	COLORREF color = RGB(colorR, colorG, colorB);
	return color;
}

//浮点数线性插值
float LerpFloat(float f0, float f1, float k)
{
	float result = f1 * k + f0 * (1 - k);
	return result;
}
//COLORREF转Vector3
Vector3 ColorToVector3(COLORREF c)
{
	Vector3 cv = Vector3(0, 0, 0);
	float f = (float)1 / 255;
	cv.x = GetRValue(c) * f;
	cv.y = GetGValue(c) * f;
	cv.z = GetBValue(c) * f;
	return cv;
}
//Vector3转COLORREF
COLORREF Vector3ToColor(Vector3 v)
{
	COLORREF c = RGB(min(v.x, 1.0f) * 255, min(v.y, 1.0f) * 255, min(v.z, 1.0f) * 255);
	return c;
}

#pragma endregion

#pragma region Math

const float PI = 3.14159265358f;
const float D2R(float deg)
{
	return deg * PI / 180;
}

//Matrix44 成员

Matrix44::Matrix44()
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			m[i][j] = 0;
		}
	}
}
Matrix44 Matrix44::Identity()
{
	Matrix44 i;
	i.m[0][0] = i.m[1][1] = i.m[2][2] = i.m[3][3] = 1;
	return i;
}


Matrix44 Matrix44::operator*(const Matrix44 &v)
{
	Matrix44 r = Matrix44();

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			r.m[i][j] = m[i][0] * v.m[0][j] + m[i][1] * v.m[1][j] + m[i][2] * v.m[2][j] + m[i][3] * v.m[3][j];
		}
	}

	return r;
}
Vector4 Matrix44::operator*(const Vector4 &v)
{
	float rx = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] * v.w;
	float ry = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] * v.w;
	float rz = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] * v.w;
	float rw = m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3] * v.w;
	return Vector4(rx, ry, rz, rw);
}

Matrix44 Matrix44::Transpose()
{
	Matrix44 mt;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mt.m[i][j] = m[j][i];
		}
	}
	return mt;
}


Matrix44 Matrix44::InverseDiagonal()
{
	Matrix44 mt = *this;
	for (int i = 0; i < 4; i++)
	{
		mt.m[i][i] = mt.m[i][i] == 0 ? 0 : 1 / mt.m[i][i];
	}
	return mt;
}

//Vector4 成员函数

Vector4::Vector4(float x = 0, float y = 0, float z = 0, float w = 0) : x(x), y(y), z(z), w(w) {}

Vector4::Vector4(Vector3 v, float w) : x(v.x), y(v.y), z(v.z), w(w) {}

Vector4 Vector4::operator+(const Vector4 &v)
{
	return Vector4(x + v.x, y + v.y, z + v.z, w + v.w);
}
Vector4 Vector4::operator-(const Vector4 &v)
{
	return Vector4(x - v.x, y - v.y, z - v.z, w + v.w);
}
float Vector4::operator*(const Vector4 &v)
{
	return x * v.x + y * v.y + z * v.z + w * v.w;
}
Vector4 Vector4::operator*(const Matrix44 &v)
{
	float rx = x * v.m[0][0] + y * v.m[1][0] + z * v.m[2][0] + w * v.m[3][0];
	float ry = x * v.m[0][1] + y * v.m[1][1] + z * v.m[2][1] + w * v.m[3][1];
	float rz = x * v.m[0][2] + y * v.m[1][2] + z * v.m[2][2] + w * v.m[3][2];
	float rw = x * v.m[0][3] + y * v.m[1][3] + z * v.m[2][3] + w * v.m[3][3];
	return Vector4(rx, ry, rz, rw);
}
Vector4 Vector4::operator*(const float &f)
{
	return Vector4(x * f, y * f, z * f, w * f);
}
Vector4 Vector4::operator^(const Vector4 &v)
{
	return Vector4(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x, 1);
}
float Vector4::GetLength()
{
	return sqrt(x * x + y * y + z * z + w * w);
}
Vector4 Vector4::Normalize()
{
	float length = GetLength();
	float factor = length == 0 ? 0 : 1 / length;
	return Vector4(x, y, z, w) * factor;
}



//Vector3 成员函数

Vector3::Vector3(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}

Vector3::Vector3(Vector4 v) : x(v.x), y(v.y), z(v.z) {}

Vector3 Vector3::operator+(const Vector3 &v)
{
	return Vector3(x + v.x, y + v.y, z + v.z);
}
Vector3 Vector3::operator-(const Vector3 &v)
{
	return Vector3(x - v.x, y - v.y, z - v.z);
}
Vector3 Vector3::operator-()
{
	return Vector3(-x, -y, -z);
}
float Vector3::operator*(const Vector3 &v)
{
	return x * v.x + y * v.y + z * v.z;
}
Vector3 Vector3::operator*(const float &f)
{
	return Vector3(x * f, y * f, z * f);
}
Vector3 Vector3::operator^(const Vector3 &v)
{
	return Vector3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
}
float Vector3::GetLength()
{
	return sqrt(x * x + y * y + z * z);
}
Vector3 Vector3::Normalize()
{
	float length = GetLength();
	float factor = length == 0 ? 0 : 1 / length;
	return Vector3(x, y, z) * factor;
}

Vector3 Vector3::Lerp(Vector3 v0, Vector3 v1, float k)
{
	return v0 * (1 - k) + v1 * k;
}

Vector3 Vector3::HadamardProduct(Vector3 v0, Vector3 v1)
{
	Vector3 v;
	v.x = v0.x * v1.x;
	v.y = v0.y * v1.y;
	v.z = v0.z * v1.z;
	return v;
}
Vector3 Vector3::Reflect(Vector3 v, Vector3 n)
{
	Vector3 r = n * (v * n) * 2 - v;
	return r;
}



Vector2 Vector2::operator+(const Vector2 &v)
{
	return Vector2(x + v.x, y + v.y);
}
Vector2 Vector2::operator-(const Vector2 &v)
{
	return Vector2(x - v.x, y - v.y);
}
float Vector2::operator*(const Vector2 &v)
{
	return x * v.x + y * v.y;
}
Vector2 Vector2::operator*(const float &f)
{
	return Vector2(x * f, y * f);
}
Vector2::Vector2(float x = 0, float y = 0) :x(x), y(y) {};
Vector2 Vector2::Lerp(Vector2 v0, Vector2 v1, float k)
{
	return v0 * (1 - k) + v1 * k;
}

#pragma endregion

#pragma region Model Struct



Vertex::Vertex(Vector3 pos = Vector3(), TexCoord texcoord = TexCoord(), Vector3 normal = Vector3(), COLORREF color = RGB(0, 0, 0)) :pos(pos), texcoord(texcoord), normal(normal), color(color) {};
Vertex Vertex::Lerp(Vertex v0, Vertex v1, float k) 
{
	Vertex r;

	r.pos = Vector3::Lerp(v0.pos, v1.pos, k);
	r.normal = Vector3::Lerp(v0.normal, v1.normal, k);
	r.color = LerpColor(v0.color, v1.color, k);

	TexCoord t;
	if (isPerspectiveCorrectOn)
	{
		//uv透视校正
		//rhw为View空间的 1 / z
		float z = max(1 / LerpFloat(v0.rhw, v1.rhw, k), 0);
		r.rhw = 1 / z;
		//r.pos = Vector3(LerpFloat(v0.pos.x, v1.pos.x, k), LerpFloat(v0.pos.y, v1.pos.y, k), z);
		t.x = max(z * LerpFloat(v0.texcoord.x * v0.rhw, v1.texcoord.x * v1.rhw, k), 0);
		t.y = max(z * LerpFloat(v0.texcoord.y * v0.rhw, v1.texcoord.y * v1.rhw, k), 0);
	}
	else
	{
		//不进行透视校正
		r.pos = Vector3::Lerp(v0.pos, v1.pos, k);
		t.x = LerpFloat(v0.texcoord.x, v1.texcoord.x, k);
		t.y = LerpFloat(v0.texcoord.y, v1.texcoord.y, k);
	}

	r.texcoord = t;

	return r;
}

Triangle::Triangle(Vertex v0, Vertex v1, Vertex v2)
{
	v[0] = v0;
	v[1] = v1;
	v[2] = v2;
};

#pragma endregion

#pragma region Camera

Camera::Camera(Vector3 pos = Vector3(0, 0, 0), Vector3 lookAt = Vector3(0, 0, 0), Vector3 up = Vector3(0, 0, 0)) :pos(pos), lookAt(lookAt), up(up) {}

#pragma endregion

#pragma region Model Loader



void Mesh::LoadOBJ(const char* filename)
{
	std::ifstream in(filename, std::ios::in);
	char curLine[100];

	if (!in.is_open())
	{
		std::cout << "Fail to open the file" << std::endl;
		return;
	}

	std::string s;
	std::string tag;
	std::vector<std::string> temp;
	while (!in.eof())
	{
		in.getline(curLine, 100);
		//使用stringstream分割字符串
		std::istringstream ss(curLine);
		while (ss >> s)
		{
			if (s == "v")
			{
				SetValue(tag, temp);
				tag = "v";
				temp.clear();
				continue;
			}
			if (s == "vt")
			{
				SetValue(tag, temp);
				tag = "vt";
				temp.clear();
				continue;
			}
			if (s == "vn")
			{
				SetValue(tag, temp);
				tag = "vn";
				temp.clear();
				continue;
			}
			if (s == "f")
			{
				SetValue(tag, temp);
				tag = "f";
				temp.clear();
				continue;
			}

			temp.push_back(s);

			//std::cout << s << std::endl;
		}
		SetValue(tag, temp);
		tag = "";
	}
	std::cout << "Load finished!" << std::endl;
}


void Mesh::SetValue(const std::string tag, std::vector<std::string> temp)
{
	if (temp.size() == 0 || tag == "")
	{
		return;
	}
	if (tag == "v")
	{
		v.push_back(Vector3(std::stof(temp[0]), std::stof(temp[1]), std::stof(temp[2])));
	}
	if (tag == "vt")
	{
		vt.push_back(TexCoord(std::stof(temp[0]), std::stof(temp[1])));

	}
	if (tag == "vn")
	{
		vn.push_back(Vector3(std::stof(temp[0]), std::stof(temp[1]), std::stof(temp[2])));

	}
	if (tag == "f")
	{
		ParseFace(temp[0], temp[1], temp[2]);
	}

}

void Mesh::ParseFace(std::string s0, std::string s1, std::string s2)
{
	std::vector<std::string> sv0, sv1, sv2;
	SplitString(sv0, s0, "/");
	SplitString(sv1, s1, "/");
	SplitString(sv2, s2, "/");

	COLORREF rgb[] = { RGB(255, 255, 0), RGB(255, 255, 0), RGB(255, 255, 0) };
	//颜色暂时根据顶点序号%3决定
	COLORREF c0 = rgb[(std::stoi(sv0[0]) - 1) % 3];
	COLORREF c1 = rgb[(std::stoi(sv1[0]) - 1) % 3];
	COLORREF c2 = rgb[(std::stoi(sv2[0]) - 1) % 3];


	Vertex v0 = Vertex(v[std::stoi(sv0[0]) - 1], vt[std::stoi(sv0[1]) - 1], vn[std::stoi(sv0[2]) - 1], c0);
	Vertex v1 = Vertex(v[std::stoi(sv1[0]) - 1], vt[std::stoi(sv1[1]) - 1], vn[std::stoi(sv1[2]) - 1], c1);
	Vertex v2 = Vertex(v[std::stoi(sv2[0]) - 1], vt[std::stoi(sv2[1]) - 1], vn[std::stoi(sv2[2]) - 1], c2);
	f.push_back(Triangle(v0, v1, v2));

}

void Mesh::SplitString(std::vector<std::string>& v, const std::string& s, const std::string& d)
{
	std::string::size_type pos1, pos2;
	pos2 = s.find(d);
	pos1 = 0;
	while (std::string::npos != pos2)
	{
		v.push_back(s.substr(pos1, pos2 - pos1));

		pos1 = pos2 + d.size();
		pos2 = s.find(d, pos1);
	}
	if (pos1 != s.length())
	{
		v.push_back(s.substr(pos1));
	}
}
#pragma endregion

#pragma region Coordinate Transform 

Transform::Transform()
{
	translate = Matrix44::Identity();
	rotation = Matrix44::Identity();
	scale = Matrix44::Identity();
	modelMatrix = Matrix44::Identity();
	viewMatrix = Matrix44::Identity();
	projectionMatrix = Matrix44::Identity();

}


void Transform::SetTranslate(float x, float y, float z)
{
	Matrix44 t = Matrix44::Identity();
	t.m[3][0] = x;
	t.m[3][1] = y;
	t.m[3][2] = z;
	translate = t;

}
void Transform::SetRotation(float x, float y, float z)
{
	Matrix44 rx = Matrix44::Identity();
	Matrix44 ry = Matrix44::Identity();
	Matrix44 rz = Matrix44::Identity();

	rx.m[1][1] = cosf(D2R(x));
	rx.m[2][2] = cosf(D2R(x));
	rx.m[2][1] = -sinf(D2R(x));
	rx.m[1][2] = sinf(D2R(x));

	ry.m[2][2] = cosf(D2R(y));
	ry.m[0][0] = cosf(D2R(y));
	ry.m[0][2] = -sinf(D2R(y));
	ry.m[2][0] = sinf(D2R(y));

	rz.m[0][0] = cosf(D2R(z));
	rz.m[1][1] = cosf(D2R(z));
	rz.m[1][0] = -sinf(D2R(z));
	rz.m[0][1] = sinf(D2R(z));

	//按照zxy顺序旋转
	rotation = rz * rx * ry;
}
void Transform::SetScale(float x, float y, float z)
{
	Matrix44 s = Matrix44();
	s.m[0][0] = x;
	s.m[1][1] = y;
	s.m[2][2] = z;
	scale = s;
}
void Transform::UpdateModelMatrix()
{
	modelMatrix = scale * rotation * translate;
}

void Transform::SetCamera(Vector3 pos, Vector3 lookAt, Vector3 up)
{
	camera = Camera(pos, lookAt, up);
}
void Transform::UpdateViewMatrix()
{
	Vector3 forward = (camera.lookAt - camera.pos).Normalize();
	Vector3 right = (camera.up ^ forward).Normalize();
	//使up与其他两向量正交
	Vector3 up = (forward ^ right).Normalize();

	Matrix44 t = Matrix44::Identity();
	t.m[3][0] = -camera.pos.x;
	t.m[3][1] = -camera.pos.y;
	t.m[3][2] = -camera.pos.z;

	Matrix44 r = Matrix44::Identity();
	r.m[0][0] = right.x;
	r.m[1][0] = right.y;
	r.m[2][0] = right.z;
	r.m[0][1] = up.x;
	r.m[1][1] = up.y;
	r.m[2][1] = up.z;
	r.m[0][2] = forward.x;
	r.m[1][2] = forward.y;
	r.m[2][2] = forward.z;

	viewMatrix = t * r;

}

void Transform::SetprojectionMatrix(float fovy, float w, float h, float zn, float zf)
{
	float aspect = w / h;
	float cot = 1 / tanf(fovy * 0.5f);
	float Q = zf / (zf - zn);

	Matrix44 p;
	p.m[0][0] = cot / aspect;
	p.m[1][1] = cot;
	p.m[2][2] = Q;
	p.m[3][2] = -Q * zn;
	p.m[2][3] = 1;
	projectionMatrix = p;
}
void Transform::UpdateMVP()
{
	MVPMatrix = modelMatrix * viewMatrix * projectionMatrix;
}

void Transform::UpdateNormalMatrix()
{
	//法向量变换到世界坐标系G = M的逆转置，不需平移，则为S逆*R
	normalMatrix = scale.InverseDiagonal() * rotation;
}

void Transform::SetViewportMatrix(float x, float y, float w, float h, float zmin = 0, float zmax = 1)
{
	Matrix44 vp = Matrix44::Identity();
	vp.m[0][0] = w * 0.5f;
	vp.m[3][0] = x + w * 0.5f;
	vp.m[1][1] = -h * 0.5f;
	vp.m[3][1] = y + h * 0.5f;
	vp.m[2][2] = zmax - zmin;
	vp.m[3][2] = zmin;
	viewportMatrix = vp;
}

void Transform::Update()
{
	UpdateModelMatrix();
	UpdateViewMatrix();
	UpdateMVP();
	UpdateNormalMatrix();
}
Vector4 Transform::PerspectiveDivision(Vector4 v)
{
	float rhw = v.w == 0 ? 0 : 1 / v.w;
	v.x *= rhw;
	v.y *= rhw;
	v.z *= rhw;
	//w结果应为1
	v.w *= rhw;
	return v;
}
Vertex Transform::Apply(Vertex v)
{
	Vertex result = v;
	Vector4 v4 = Vector4(v.pos, 1);
	Vector4 r = v4 * MVPMatrix;
	result.pos = Vector3(PerspectiveDivision(r) * viewportMatrix);
	result.rhw = 1 / r.w;

	result.wPos = Vector3(v4 * modelMatrix);
	result.normal = Vector3(Vector4(v.normal, 0) * normalMatrix).Normalize();

	return result;
}

#pragma endregion

#pragma region Buffer

RenderBuffer::RenderBuffer(int w, int h)
{
	width = w;
	height = h;
	zBuffer = std::vector<float>(w * h, 1);
	frameBuffer = std::vector<COLORREF>(w * h, RGB(255, 255, 255));
}

void RenderBuffer::Resize(int w, int h)
{
	width = w;
	height = h;
	zBuffer = std::vector<float>(w * h, 1);
	frameBuffer = std::vector<COLORREF>(w * h, RGB(255, 255, 255));
}



#pragma endregion

#pragma region Buffer Draw

//绘制直线-DDA算法
void RenderBuffer::LineDDA(Vertex v0, Vertex v1)
{

	float k = float(v1.pos.y - v0.pos.y) / (v1.pos.x - v0.pos.x);
	bool isSteep = false;

	if (fabs(k) > 1)
	{
		k = 1 / k;
		std::swap(v0.pos.y, v0.pos.x);
		std::swap(v1.pos.y, v1.pos.x);
		isSteep = true;
	}

	if (v0.pos.x > v1.pos.x)
	{
		std::swap(v0, v1);
	}

	float y = v0.pos.y - k;
	//当前比例
	float ratio;
	Vertex v;
	for (int x = v0.pos.x; x <= v1.pos.x; x++)
	{
		y += k;
		ratio = (x - v0.pos.x) / (v1.pos.x - v0.pos.x);
		v = Vertex::Lerp(v0, v1, ratio);
		if (isSteep)
		{
			if (v.pos.z >= 0 && v.pos.z < zBuffer[(int)roundf(y) + x * width])
			{
				frameBuffer[(int)roundf(y) + x * width] = v.color;
				zBuffer[(int)roundf(y) + x * width] = v.pos.z;
			}
		}
		else
		{
			if (v.pos.z >= 0 && v.pos.z < zBuffer[x + (int)roundf(y) * width])
			{
				frameBuffer[x + (int)roundf(y) * width] = v.color;
				zBuffer[x + (int)roundf(y) * width] = v.pos.z;
			}
		}
	}
}

//绘制直线-Bresenham算法
void RenderBuffer::LineBresenham(Vertex v0, Vertex v1)
{
	float k = float(v1.pos.y - v0.pos.y) / (v1.pos.x - v0.pos.x);
	bool isSteep = false;

	//根据斜率1分界, xy互换
	if (fabs(k) > 1)
	{
		k = 1 / k;
		std::swap(v0.pos.y, v0.pos.x);
		std::swap(v1.pos.y, v1.pos.x);
		isSteep = true;
	}

	//起止点互换，从x较小点开始
	if (v0.pos.x > v1.pos.x)
	{
		std::swap(v0, v1);
	}

	//k < 0时，计算F增量时dx取-dx
	int dx = k > 0 ? (v1.pos.x - v0.pos.x) : (v0.pos.x - v1.pos.x);
	int dy = v1.pos.y - v0.pos.y;
	int F = 2 * dy - dx;
	int y = v0.pos.y;


	float ratio;
	Vertex v;
	for (int x = v0.pos.x; x <= v1.pos.x; x++)
	{
		//k > 0时，判断F > 0; k < 0时，判断F < 0
		if (k * F > 0)
		{
			k > 0 ? y++ : y--;
			F += 2 * (dy - dx);
		}
		else
		{
			F += 2 * dy;
		}

		//根据斜率1分界, xy绘制时换回
		ratio = (x - v0.pos.x) / (v1.pos.x - v0.pos.x);
		v = Vertex::Lerp(v0, v1, ratio);
		if (isSteep)
		{
			if (v.pos.z >= 0 && v.pos.z < zBuffer[(int)roundf(y) + x * width])
			{
				frameBuffer[(int)roundf(y) + x * width] = v.color;
				zBuffer[(int)roundf(y) + x * width] = v.pos.z;
			}
		}
		else
		{
			if (v.pos.z >= 0 && v.pos.z < zBuffer[x + (int)roundf(y) * width])
			{
				frameBuffer[x + (int)roundf(y) * width] = v.color;
				zBuffer[x + (int)roundf(y) * width] = v.pos.z;
			}
		}
	}
}

void RenderBuffer::FillTriangle(Vertex v0, Vertex v1, Vertex v2)
{
	if (v0.pos.y == v1.pos.y && v0.pos.y == v2.pos.y)
	{
		std::cout << "y1 = y2 = y3" << std::endl;
		return;
	}
	//根据y对顶点冒泡排序v0, v1, v2
	if (v0.pos.y > v1.pos.y)
	{
		std::swap(v0, v1);
	}
	if (v1.pos.y > v2.pos.y)
	{
		std::swap(v1, v2);
	}
	if (v0.pos.y > v1.pos.y)
	{
		std::swap(v0, v1);
	}

	//z和color浮点数插值
	float step01 = (v1.pos.y - v0.pos.y) == 0 ? 0 : 1 / (v1.pos.y - v0.pos.y);
	float step02 = (v2.pos.y - v0.pos.y) == 0 ? 0 : 1 / (v2.pos.y - v0.pos.y);
	float step12 = (v2.pos.y - v1.pos.y) == 0 ? 0 : 1 / (v2.pos.y - v1.pos.y);

	float k01 = 0;
	float k02 = 0;
	float k12 = 0;

	//Vertex v01, v02, v03;

	int x01 = (int)v0.pos.x;
	int x02 = (int)v0.pos.x;
	int x12 = (int)v1.pos.x;

	//Scanline-Zbuffer算法
	int scanlineY = (int)v0.pos.y;
	for (; scanlineY < v1.pos.y; scanlineY++)
	{
		if (isTex)
		{
			ScanLineTex(scanlineY, x01, x02, Vertex::Lerp(v0, v1, k01), Vertex::Lerp(v0, v2, k02));
		}
		else
		{
			ScanLine(scanlineY, x01, x02, Vertex::Lerp(v0, v1, k01), Vertex::Lerp(v0, v2, k02));
		}

		k01 += step01;
		k01 = min(k01, 1.0f);
		k02 += step02;
		k02 = min(k02, 1.0f);

		//k01 = (scanlineY - v0.pos.y) / (v1.pos.y - v0.pos.y);
		//k02 = (scanlineY - v0.pos.y) / (v2.pos.y - v0.pos.y);

		x01 = (int)LerpFloat(v0.pos.x, v1.pos.x, k01);
		x02 = (int)LerpFloat(v0.pos.x, v2.pos.x, k02);
	}
	for (; scanlineY <= v2.pos.y; scanlineY++)
	{
		if (isTex)
		{
			ScanLineTex(scanlineY, x12, x02, Vertex::Lerp(v1, v2, k12), Vertex::Lerp(v0, v2, k02));
		}
		else
		{
			ScanLine(scanlineY, x12, x02, Vertex::Lerp(v1, v2, k12), Vertex::Lerp(v0, v2, k02));
		}

		k12 += step12;
		k12 = min(k12, 1.0f);
		k02 += step02;
		k02 = min(k02, 1.0f);

		//k12 = (scanlineY - v1.pos.y) / (v2.pos.y - v1.pos.y);
		//k02 = (scanlineY - v0.pos.y) / (v2.pos.y - v0.pos.y);
		
		x12 = (int)LerpFloat(v1.pos.x, v2.pos.x, k12);
		x02 = (int)LerpFloat(v0.pos.x, v2.pos.x, k02);
	}
}

void RenderBuffer::ScanLine(int y, int x0, int x1, Vertex v0, Vertex v1)
{
	float k;
	if (x0 > x1)
	{
		std::swap(x0, x1);
		std::swap(v0, v1);
	}


	COLORREF c;
	Vertex vertex;
	int u, v;
	for (int i = x0; i <= x1; i++)
	{
		k = x1 == x0 ? 0 : (float)(i - x0) / (x1 - x0);
		vertex = Vertex::Lerp(v0, v1, k);
		c = Vector3ToColor(light.Apply(vertex, trans.camera) + ambientColor * ambientStrength);
		FillPixel(i, y, LerpFloat(v0.pos.z, v1.pos.z, k), c);
	}
}

void RenderBuffer::ScanLineTex(int y, int x0, int x1, Vertex v0, Vertex v1)
{
	float k;
	if (x0 > x1)
	{
		std::swap(x0, x1);
		std::swap(v0, v1);
	}
	COLORREF c;
	Vertex vertex;
	int u, v;
	for (int i = x0; i <= x1; i++)
	{
		k = x1 == x0 ? 0 : (float)(i - x0) / (x1 - x0);
		vertex = Vertex::Lerp(v0, v1, k);
		u = min((int)floorl(vertex.texcoord.x * texWidth), texWidth - 1);
		v = min((int)floorl(vertex.texcoord.y * texHeight), texHeight - 1);
		vertex.color = tex.GetPixel(u, v);
		c = Vector3ToColor(light.Apply(vertex, trans.camera) + ambientColor * ambientStrength);
		FillPixel(i, y, LerpFloat(v0.pos.z, v1.pos.z, k), c);
	}
}

void RenderBuffer::FillPixel(int x, int y, float z, COLORREF c)
{
	if (y >= height || y < 0 || x >= width || x < 0)
	{
		return;
	}
	if (z >= 0 && z < zBuffer[x + y * width])
	{
		zBuffer[x + y * width] = z;
		frameBuffer[x + y * width] = c;
	}
}

void RenderBuffer::DrawFrame(HDC hdc)
{
	/*for (int i = 0; i < width * height; i++)
	{
		DrawPixel(hdc, i % width, i / width, frameBuffer[i]);
	}*/
	BITMAPINFO bi;
	bi.bmiHeader.biSize = sizeof(BITMAPINFO);
	bi.bmiHeader.biWidth = width;
	//从上向下扫描需要设为-height
	bi.bmiHeader.biHeight = -height;
	bi.bmiHeader.biPlanes = 1;
	bi.bmiHeader.biBitCount = sizeof(COLORREF) * 8;
	bi.bmiHeader.biCompression = BI_RGB;
	bi.bmiHeader.biSizeImage = (DWORD)(width * height * sizeof(COLORREF));
	bi.bmiHeader.biXPelsPerMeter = bi.bmiHeader.biYPelsPerMeter = 72;
	bi.bmiHeader.biClrUsed = 0;
	bi.bmiHeader.biClrImportant = 0;


	auto error_code = SetDIBitsToDevice(hdc, 0, 0, width, height, 0, 0, 0, height, &frameBuffer[0], &bi, DIB_RGB_COLORS);
};


TexCoord RenderBuffer::LerpTexCoord(Vertex v0, Vertex v1, float k)
{
	float z = 1 / (1 / v0.pos.z + k * (1 / v1.pos.z - 1 / v0.pos.z));
	TexCoord t;
	t.x = (v0.texcoord.x / v0.pos.z + k * (v1.texcoord.x / v1.pos.z - v0.texcoord.x / v0.pos.z)) * z;
	t.y = (v0.texcoord.y / v0.pos.z + k * (v1.texcoord.y / v1.pos.z - v0.texcoord.y / v0.pos.z)) * z;
	return t;
}

#pragma endregion

#pragma region Draw
//
//
//
//绘制单个像素点
void DrawPixel(HDC hdc, int x, int y, COLORREF color)
{
	SetPixel(hdc, x, y, color);
}
//
////绘制直线-DDA算法
//void DrawLineDDA(HDC hdc, int x0, int y0, int x1, int y1, COLORREF color0, COLORREF color1)
//{
//
//	float k = float(y1 - y0) / (x1 - x0);
//	bool isSteep = false;
//
//	if (fabs(k) > 1)
//	{
//		k = 1 / k;
//		std::swap(x0, y0);
//		std::swap(x1, y1);
//		isSteep = true;
//	}
//
//	if (x0 > x1)
//	{
//		std::swap(x0, x1);
//		std::swap(y0, y1);
//		std::swap(color0, color1);
//	}
//
//	float y = y0 - k;
//	for (int x = x0; x <= x1; x++)
//	{
//		y += k;
//		if (isSteep)
//		{
//			DrawPixel(hdc, (int)roundf(y), x, LerpColor(color0, color1, float(x - x0) / (x1 - x0)));
//		}
//		else
//		{
//			DrawPixel(hdc, x, (int)roundf(y), LerpColor(color0, color1, float(x - x0) / (x1 - x0)));
//		}
//	}
//}
//
//绘制直线-Bresenham算法
void DrawLineBresenham(HDC hdc, int x0, int y0, int x1, int y1, COLORREF color0, COLORREF color1)
{
	float k = float(y1 - y0) / (x1 - x0);
	bool isSteep = false;

	//根据斜率1分界, xy互换
	if (fabs(k) > 1)
	{
		k = 1 / k;
		std::swap(x0, y0);
		std::swap(x1, y1);
		isSteep = true;
	}

	//起止点互换，从x较小点开始
	if (x0 > x1)
	{
		std::swap(x0, x1);
		std::swap(y0, y1);
		std::swap(color0, color1);
	}

	//k < 0时，计算F增量时dx取-dx
	int dx = k > 0 ? (x1 - x0) : (x0 - x1);
	int dy = y1 - y0;
	int F = 2 * dy - dx;
	int y = y0;

	for (int x = x0; x <= x1; x++)
	{
		//k > 0时，判断F > 0; k < 0时，判断F < 0
		if (k * F > 0)
		{
			k > 0 ? y++ : y--;
			F += 2 * (dy - dx);
		}
		else
		{
			F += 2 * dy;
		}

		//根据斜率1分界, xy绘制时换回
		if (isSteep)
		{
			DrawPixel(hdc, y, x, LerpColor(color0, color1, float(x - x0) / (x1 - x0)));
		}
		else
		{
			DrawPixel(hdc, x, y, LerpColor(color0, color1, float(x - x0) / (x1 - x0)));
		}
	}
}

////三角形填充
//void FillBottomFlatTriangle(HDC hdc, Vertex v0, Vertex v1, Vertex v2)
//{
//	if (v1.pos.y != v2.pos.y)
//	{
//		std::cout << "not a bottom flat triangle" << std::endl;
//		return;
//	}
//
//	float k01 = (v1.pos.y - v0.pos.y) / (v1.pos.x - v0.pos.x);
//	float k02 = (v2.pos.y - v0.pos.y) / (v2.pos.x - v0.pos.x);
//	float rk01 = 1 / k01;
//	float rk02 = 1 / k02;
//
//	float x1 = v0.pos.x;
//	float x2 = v0.pos.x;
//
//	for (int scanlineY = v0.pos.y; scanlineY <= v2.pos.y; scanlineY++)
//	{
//		DrawLineBresenham(hdc, (int)x1, scanlineY, (int)x2, scanlineY, RGB(255, 0, 0), RGB(0, 255, 255));
//		x1 += rk01;
//		x2 += rk02;
//	}
//}
//
//void FillTopFlatTriangle(HDC hdc, Vertex v0, Vertex v1, Vertex v2)
//{
//	if (v1.pos.y != v2.pos.y)
//	{
//		std::cout << "not a top flat triangle" << std::endl;
//	}
//
//	float k01 = (v1.pos.y - v0.pos.y) / (v1.pos.x - v0.pos.x);
//	float k02 = (v2.pos.y - v0.pos.y) / (v2.pos.x - v0.pos.x);
//	float rk01 = 1 / k01;
//	float rk02 = 1 / k02;
//
//	float x1 = v0.pos.x;
//	float x2 = v0.pos.x;
//
//	for (int scanlineY = v0.pos.y; scanlineY > v2.pos.y; scanlineY--)
//	{
//		DrawLineBresenham(hdc, (int)x1, scanlineY, (int)x2, scanlineY, RGB(255, 0, 0), RGB(0, 255, 255));
//		x1 -= rk01;
//		x2 -= rk02;
//	}
//}
//
//void FillTriangle(HDC hdc, Vertex v0, Vertex v1, Vertex v2)
//{
//	if (v0.pos.y == v1.pos.y && v0.pos.y == v2.pos.y)
//	{
//		std::cout << "y1 = y2 = y3" << std::endl;
//		return;
//	}
//	//根据y对顶点冒泡排序
//	if (v0.pos.y > v1.pos.y)
//	{
//		std::swap(v0, v1);
//	}
//	if (v1.pos.y > v2.pos.y)
//	{
//		std::swap(v1, v2);
//	}
//	if (v0.pos.y > v1.pos.y)
//	{
//		std::swap(v0, v1);
//	}
//
//	//创建中间分割点v3
//	float k = (v1.pos.y - v0.pos.y) / (v2.pos.y - v0.pos.y);
//
//	Vertex v3 = Vertex::Lerp(v0, v2, k);
//	FillBottomFlatTriangle(hdc, v0, v1, v3);
//	FillTopFlatTriangle(hdc, v2, v1, v3);
//
//}
//
#pragma endregion

#pragma region Light


DirectionalLight::DirectionalLight() {};

Vector3 DirectionalLight::Apply(Vertex v, Camera camera)
{
	direction = direction.Normalize();
	float diffuseFactor = max(-direction * v.normal, 0);
	Vector3 diffuse = Vector3::HadamardProduct(ColorToVector3(v.color), color) * diffuseFactor;

	Vector3 reflect = Vector3::Reflect(-direction, v.normal).Normalize();
	float specularFactor = pow(max(reflect * (camera.pos - v.wPos).Normalize(), 0), shininess);
	Vector3 specular = Vector3::HadamardProduct(ColorToVector3(v.color), color) * specularFactor;

	Vector3 result = diffuse + specular;
	return result;
}


#pragma endregion



void Init(HWND hWnd)
{
	//获取窗口宽高
	RECT clientRect;
	GetClientRect(hWnd, &clientRect);
	int width = clientRect.right - clientRect.left;
	int height = clientRect.bottom - clientRect.top;

	//新建缓冲区
	rb = RenderBuffer(width, height);


	//读取模型
	mesh.LoadOBJ(filename);

	//设置变换
	trans.SetRotation(rotationX, rotationY, 0);
	trans.SetCamera(Vector3(cameraX, 0, cameraZ), Vector3(cameraLookAtX, 0, 0), Vector3(0, 1, 0));
	trans.SetprojectionMatrix(PI * 0.5f, width, height, 1, 1000);
	trans.SetViewportMatrix(0, 0, width, height);
	trans.Update();

	//设置方向光
	light.direction = lightDirection;
	light.color = lightColor;

	
	HRESULT h = tex.Load(texFilename);

	texWidth = tex.GetWidth();
	texHeight = tex.GetHeight();

	std::cout << "Init" << std::endl;
	std::cout << "Tex Width:" << texWidth << " Tex Height:" << texHeight << std::endl;
}

void Render(HWND hWnd)
{
	PAINTSTRUCT ps;
	HDC hdc = BeginPaint(hWnd, &ps);

	//获取窗口宽高
	RECT clientRect;
	GetClientRect(hWnd, &clientRect);
	int width = clientRect.right - clientRect.left;
	int height = clientRect.bottom - clientRect.top;

	rb.Resize(width, height);


	////读取模型
	//Mesh mesh;
	//mesh.LoadOBJ(filename);

	////设置变换
	//Transform trans;
	//trans.SetRotation(30, 30, 0);
	//trans.SetCamera(Vector3(0, 0, 3), Vector3(0, 0, 0), Vector3(0, 1, 0));
	//trans.SetprojectionMatrix(PI * 0.5f, width, height, 1, 1000);
	//trans.SetViewportMatrix(0, 0, width, height);

	trans.Update();

	
	Vertex v0, v1, v2;

	for (unsigned int i = 0; i < mesh.f.size(); i++)
	{
		Triangle t = mesh.f[i];
		/*Vector3 v0 = trans.Apply(t.v[0].pos);
		Vector3 v1 = trans.Apply(t.v[1].pos);
		Vector3 v2 = trans.Apply(t.v[2].pos);

		t.v[0].pos = v0;
		t.v[1].pos = v1;
		t.v[2].pos = v2;*/

		v0 = trans.Apply(t.v[0]);
		v1 = trans.Apply(t.v[1]);
		v2 = trans.Apply(t.v[2]);

		//v0.color = Vector3ToColor(light.Apply(v0, trans.camera) + ambientColor * ambientStrength);
		//v1.color = Vector3ToColor(light.Apply(v1, trans.camera) + ambientColor * ambientStrength);
		//v2.color = Vector3ToColor(light.Apply(v2, trans.camera) + ambientColor * ambientStrength);

		/*DrawLineBresenham(hdc, (int)v0.pos.x, (int)v0.pos.y, (int)v1.pos.x, (int)v1.pos.y, RGB(255, 0, 0), RGB(0, 255, 0));
		DrawLineBresenham(hdc, (int)v1.pos.x, (int)v1.pos.y, (int)v2.pos.x, (int)v2.pos.y, RGB(0, 255, 0), RGB(0, 0, 255));
		DrawLineBresenham(hdc, (int)v2.pos.x, (int)v2.pos.y, (int)v0.pos.x, (int)v0.pos.y, RGB(0, 0, 255), RGB(255, 0, 0));*/

		rb.FillTriangle(v0, v1, v2);
	}

	rb.DrawFrame(hdc);
	EndPaint(hWnd, &ps);

	std::cout << "camera pos: (" << trans.camera.pos.x << "," << trans.camera.pos.y << "," << trans.camera.pos.z << ")" << std::endl;
	std::cout << "rotation: (" << rotationX << "," << rotationY << ")" << std::endl;

}


void GameLoop(HWND hWnd)
{
	//std::cout << "game loop" << std::endl;

}


#pragma region WinMain


int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
	_In_opt_ HINSTANCE hPrevInstance,
	_In_ LPWSTR    lpCmdLine,
	_In_ int       nCmdShow)
{
	UNREFERENCED_PARAMETER(hPrevInstance);
	UNREFERENCED_PARAMETER(lpCmdLine);

	// TODO: 在此处放置代码。
	AllocConsole();
	freopen("CONOUT$", "w", stdout);
	std::cout << "Console Open" << std::endl;



	// 初始化全局字符串
	LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
	LoadStringW(hInstance, IDC_SOFTRENDERER, szWindowClass, MAX_LOADSTRING);
	MyRegisterClass(hInstance);

	// 执行应用程序初始化:
	if (!InitInstance(hInstance, nCmdShow))
	{
		return FALSE;
	}

	HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_SOFTRENDERER));

	MSG msg;

	// 主消息循环:
	while (GetMessage(&msg, nullptr, 0, 0))
	{
		if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
		GameLoop(msg.hwnd);
	}

	return (int)msg.wParam;
}

//
//  函数: MyRegisterClass()
//
//  目标: 注册窗口类。
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
	WNDCLASSEXW wcex;

	wcex.cbSize = sizeof(WNDCLASSEX);

	wcex.style = CS_HREDRAW | CS_VREDRAW;
	wcex.lpfnWndProc = WndProc;
	wcex.cbClsExtra = 0;
	wcex.cbWndExtra = 0;
	wcex.hInstance = hInstance;
	wcex.hIcon = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_SOFTRENDERER));
	wcex.hCursor = LoadCursor(nullptr, IDC_ARROW);
	wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
	wcex.lpszMenuName = MAKEINTRESOURCEW(IDC_SOFTRENDERER);
	wcex.lpszClassName = szWindowClass;
	wcex.hIconSm = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

	return RegisterClassExW(&wcex);
}

//
//   函数: InitInstance(HINSTANCE, int)
//
//   目标: 保存实例句柄并创建主窗口
//
//   注释:
//
//        在此函数中，我们在全局变量中保存实例句柄并
//        创建和显示主程序窗口。
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
	hInst = hInstance; // 将实例句柄存储在全局变量中

	HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
		CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

	//初始化
	Init(hWnd);

	if (!hWnd)
	{
		return FALSE;
	}

	ShowWindow(hWnd, nCmdShow);
	UpdateWindow(hWnd);


	return TRUE;
}




//
//  函数: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  目标: 处理主窗口的消息。
//
//  WM_COMMAND  - 处理应用程序菜单
//  WM_PAINT    - 绘制主窗口
//  WM_DESTROY  - 发送退出消息并返回
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	switch (message)
	{
	case WM_COMMAND:
	{
		int wmId = LOWORD(wParam);
		// 分析菜单选择:
		switch (wmId)
		{
		case IDM_ABOUT:
			DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
			break;
		case IDM_EXIT:
			DestroyWindow(hWnd);
			break;
		default:
			return DefWindowProc(hWnd, message, wParam, lParam);
		}
	}
	break;

	case WM_KEYDOWN:
		switch (wParam)
		{
		case VK_UP:
			rotationX -= rotationStep;
			trans.SetRotation(rotationX, rotationY, 0);
			InvalidateRect(hWnd, NULL, false);
			std::cout << "XR-" << std::endl;
			break;
		case VK_DOWN:
			rotationX += rotationStep;
			trans.SetRotation(rotationX, rotationY, 0);
			InvalidateRect(hWnd, NULL, false);
			std::cout << "XR+" << std::endl;
			break;
		case VK_LEFT:
			rotationY += rotationStep;
			trans.SetRotation(rotationX, rotationY, 0);
			InvalidateRect(hWnd, NULL, false);
			std::cout << "YR+" << std::endl;
			break;
		case VK_RIGHT:
			rotationY -= rotationStep;
			trans.SetRotation(rotationX, rotationY, 0);
			InvalidateRect(hWnd, NULL, false);
			std::cout << "YR-" << std::endl;
			break;
		case 0x57:
			cameraZ -= cameraStep;
			trans.SetCamera(Vector3(cameraX, 0, cameraZ), Vector3(cameraLookAtX, 0, 0), Vector3(0, 1, 0));
			InvalidateRect(hWnd, NULL, false);
			std::cout << "W" << std::endl;
			break;
		case 0x41:
			cameraX -= cameraStep;
			trans.SetCamera(Vector3(cameraX, 0, cameraZ), Vector3(cameraLookAtX, 0, 0), Vector3(0, 1, 0));
			InvalidateRect(hWnd, NULL, false);
			std::cout << "A" << std::endl;
			break;
		case 0x53:
			cameraZ += cameraStep;
			trans.SetCamera(Vector3(cameraX, 0, cameraZ), Vector3(cameraLookAtX, 0, 0), Vector3(0, 1, 0));
			InvalidateRect(hWnd, NULL, false);
			std::cout << "S" << std::endl;
			break;
		case 0x44:
			cameraX += cameraStep;
			trans.SetCamera(Vector3(cameraX, 0, cameraZ), Vector3(cameraLookAtX, 0, 0), Vector3(0, 1, 0));
			InvalidateRect(hWnd, NULL, false);
			std::cout << "D" << std::endl;
			break;
		case 0x51:
			cameraLookAtX += cameraLookAtStep;
			trans.SetCamera(Vector3(cameraX, 0, cameraZ), Vector3(cameraLookAtX, 0, 0), Vector3(0, 1, 0));
			InvalidateRect(hWnd, NULL, false);
			std::cout << "Q" << std::endl;
			break;
		case 0x45:
			cameraLookAtX -= cameraLookAtStep;
			trans.SetCamera(Vector3(cameraX, 0, cameraZ), Vector3(cameraLookAtX, 0, 0), Vector3(0, 1, 0));
			InvalidateRect(hWnd, NULL, false);
			std::cout << "E" << std::endl;
			break;
		}
		return 0;
	case WM_PAINT:
	{
		//绘制函数
		Render(hWnd);
		std::cout << "Paint" << std::endl;
	}
	break;
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
	return 0;
}




// “关于”框的消息处理程序。
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
	UNREFERENCED_PARAMETER(lParam);
	switch (message)
	{
	case WM_INITDIALOG:
		return (INT_PTR)TRUE;

	case WM_COMMAND:
		if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
		{
			EndDialog(hDlg, LOWORD(wParam));
			return (INT_PTR)TRUE;
		}
		break;
	}
	return (INT_PTR)FALSE;
}


#pragma endregion