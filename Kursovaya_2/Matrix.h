#pragma once
#include <set>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

//Структура хранения координат точек сетки
struct Point
{
	Point(double x, double y)
	{
		this->x = x;
		this->y = y;
	}

	Point()
	{
		x = 0;
		y = 0;
	}

	double x, y;

	Point operator+(Point& a)
	{
		return Point(a.x + x, a.y + y);
	}

	Point operator*(double constant)
	{
		return Point(x * constant, y * constant);
	}

	Point operator/(double constant)
	{
		return Point(x / constant, y / constant);
	}
};



//Структура хранящая информацию о треугольнике
struct Triangle
{
	Triangle(vector<int> Index_first, vector<int> Index_second)
	{

		this->Index_first = Index_first;
		this->Index_second = Index_second;
	}
	vector<int> Index_first, Index_second;
};

//Структура хранения глобальной матрицы
struct GlobalMatrix
{
	vector<double> ggl, ggu, di;
	vector<int> ig, jg;
};

struct Edge
{
	//int v1, v2, v3, v4;
	vector<int> v;
	int valueNo;
};

//Структура хранящая значения функций альфа
struct Alfa
{
	Alfa(double alfa1, double alfa2)
	{
		this->alfa1 = alfa1;
		this->alfa2 = alfa2;
	}

	Alfa()
	{
		alfa1 = 0;
		alfa2 = 0;
	}

	double alfa1, alfa2;
};

//Структура хранящая частные производные шаблонных функций
struct GradTemp
{
	int gradNo;
	double coeff;
	int v1, v2, v3;
};

//Структура данных произведения градиентов
struct GradMult
{
	GradMult(int grad1, int grad2, double coeff)
	{
		this->grad1 = grad1;
		this->grad2 = grad2;
		this->coeff = coeff;
	}
	int grad1, grad2;
	double coeff;
};

//Структура хранящая шаблонные функции
struct PsiTemp
{
	double coeff;
	int v1, v2, v3;
};

class Matrix
{
public:

	Matrix(int node_count, int triangle_count, vector<vector<GradTemp>> G_temp, vector<vector<PsiTemp>> M_temp)
	{
		this->node_count = node_count; //number of units
		this->triangle_count = triangle_count; //number of triangles
		this->G_temp = G_temp;
		this->M_temp = M_temp;

		memory_allocation();
		read_matrix();
		portrait();
		gen_global_G();
		gen_global_M();
	}

	vector<Point> get_units() { return Unit; }
	GlobalMatrix get_matrix_M() { return M_glob; }
	GlobalMatrix get_matrix_G() { return G_glob; }

	vector<double> get_global_vector() { return F_global; }

	void gen_global_vector(double t_now);

private:

	

	vector<Point> Unit;						// координаты вершин
	vector<Alfa> Alfa_vec;					// коэффициенты фунции разложения
	vector<double> F_global;				// вектор глобальной правой части
	vector<double> F_local;					// вектор локальной правой части
	vector<double> arg_lambda;				// коэффициенты проводимости
	vector<vector<GradTemp>> G_temp;
	vector<vector<PsiTemp>> M_temp;
	vector<vector<double>> M_coeff, G_local;	// локальная матрица массы
	vector<vector<vector<GradMult>>> G_coeff; // локальная матрица жесткости
	vector<Triangle> Ind_tri;			// номера элементов треугольника
	GlobalMatrix M_glob, G_glob;		// глобальные матрицы

	int node_count, triangle_count;

	//Функции формирования глобальной матрицы
	void gen_M();
	void gen_G();
	void gen_global_G();
	void gen_global_M();
	void gen_local_G(Triangle& tri, double lambda);
	void gen_local_vector(Triangle& tri, double t_now);
	

	double fun(double x, double y, double t);
	double factorial(int v);
	double detD(vector<Point>& P);
	
	
	void portrait();
	void read_matrix();
	void gen_Alfa(vector<Point>& P, double det);
	void add_to_glob(Triangle& tri, GlobalMatrix& Glob_mat, vector<vector<double>> Local_mat);
	void memory_allocation();
};

