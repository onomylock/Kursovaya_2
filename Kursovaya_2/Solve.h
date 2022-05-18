#pragma once
#include "Matrix.h"
#include <set>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

//Структура хранящая коэффициенты для каждого элемента
struct Argument
{
	Argument(double sigma, double hi)
	{
		this->sigma = sigma;
		this->hi = hi;
	}

	Argument()
	{
		sigma = 1;
		hi = 1;
	}

	double sigma, hi;

};

class Solve
{
public:

	Solve(int node_count, int triangle_count, int time_count, double t_start, double t_end, double t_h)
	{
		this->node_count = node_count;
		this->triangle_count = triangle_count;
		this->time_count = time_count;
		this->t_start = t_start;
		this->t_end = t_end;
		this->t_h = t_h;
		read_data();
	};
	void method();

	vector<vector<double>> get_Q() { return Q; }
	vector<vector<double>> get_U();

private:

	vector<double> Edge_basis = { 0.125, 0.375, 0.375, 0.125 };

	int node_count, triangle_count, time_count, boun_num1 = 0, boun_num2 = 0;
	double t_start, t_end, t_h;
	vector<double> F;
	vector<Argument> Arg;
	GlobalMatrix Global, G, M;
	vector<vector<double>> Q, U;
	vector<Edge> S1;
	vector<Edge> S2;
	vector<Point> Unit;

	vector<vector<GradTemp>> G_temp = {
	{ { 0, 13.5, 2, 0, 0 }, { 0, -9, 1, 0, 0 }, { 0, 1, 0, 0, 0 } }, //1
	{ { 1, 13.5, 0, 2, 0 }, { 1, -9, 0, 1, 0 }, { 1, 1, 0, 0, 0 } }, //2
	{ { 2, 13.5, 0, 0, 2 }, { 2, -9, 0, 0, 1 }, { 2, 1, 0, 0, 0 } }, //3
	{ { 0, 27, 1, 1, 0 }, { 1, 13.5, 2, 0, 0 }, { 1, -4.5, 1, 0, 0 }, { 0, -4.5, 0, 1, 0 } },//4
	{ { 1, 27, 1, 1, 0 }, { 0, 13.5, 0, 2, 0 }, { 1, -4.5, 1, 0, 0 }, { 0, -4.5, 0, 1, 0 } },//5	
	{ { 1, 27, 0, 1, 1 }, { 2, 13.5, 0, 2, 0 }, { 2, -4.5, 0, 1, 0 }, { 1, -4.5, 0, 0, 1 } },//6
	{ { 2, 27, 0, 1, 1 }, { 1, 13.5, 0, 0, 2 }, { 2, -4.5, 0, 1, 0 }, { 1, -4.5, 0, 0, 1 } },//7
	{ { 2, 27, 1, 0, 1 }, { 0, 13.5, 0, 0, 2 }, { 2, -4.5, 1, 0, 0 }, { 0, -4.5, 0, 0, 1 } },//8
	{ { 0, 27, 1, 0, 1 }, { 2, 13.5, 2, 0, 0 }, { 2, -4.5, 1, 0, 0 }, { 0, -4.5, 0, 0, 1 } },//9
	{ { 0, 27, 0, 1, 1 }, { 1, 27, 1, 0, 1 }, { 2, 27, 1, 1, 0 } } //10
	};

	vector<vector<PsiTemp>> M_temp = {
	{ { 4.5, 3, 0, 0 }, { -4.5, 2, 0, 0 }, { 1, 1, 0, 0 } }, //1
	{ { 4.5, 0, 3, 0 }, { -4.5, 0, 2, 0 }, { 1, 0, 1, 0 } }, //2
	{ { 4.5, 0, 0, 3 }, { -4.5, 0, 0, 2 }, { 1, 0, 0, 1 } }, //3
	{ { 13.5, 2, 1, 0 }, { -4.5, 1, 1, 0 } }, //4
	{ { 13.5, 1, 2, 0 }, { -4.5, 1, 1, 0 } }, //5
	{ { 13.5, 0, 2, 1 }, { -4.5, 0, 1, 1 } }, //6
	{ { 13.5, 0, 1, 2 }, { -4.5, 0, 1, 1 } }, //7
	{ { 13.5, 1, 0, 2 }, { -4.5, 1, 0, 1 } }, //8
	{ { 13.5, 2, 0, 1 }, { -4.5, 1, 0, 1 } }, //9
	{ { 27, 1, 1, 1 } } //10
	};

	void allocate_memory();
	void bound1(double t_now);
	void bound2();
	void read_data();
	void gen_matrix(double h, int iter);
	void gen_vector(vector<double> B, vector<double> B_last, double h, int iter);

	double theta_value(int num);
	double U_true(double x, double y, double t);
	double S1_gen(double x, double y, double t, int bound1);
	double distance(Point a, Point b);

	//Функции решения СЛАУ

	void LOS(vector<double>& Q);
	void LUsq(vector<double>& ggl, vector<double>& ggu, vector<double>& di);
	void mult_vec_matrix(GlobalMatrix Mat, vector<double> param, vector<double>& res);
	void mult_vec_matrix(vector<double> param, vector<double>& res);
	void forward(vector<double> ggl, vector<double> di, vector<double> F, vector<double>& X);
	void reverse(vector<double> ggu, vector<double> di, vector<double> F, vector<double>& X);
	double scalar_product(vector<double> A1, vector<double> A2);
};
