#include "Solve.h"
#include "Matrix.h"

double Solve::U_true(double x, double y, double t)
{
	//return pow(x, 4);
	//return pow(x, 3);
	//return x;
	return pow(x, 2);
	//return 3;
	//return pow(x, 4) * pow(t, 3);
	//return sin(x);
	//return pow(t, 2) * pow(x, 3);
	//return pow(x, 4);
}

double Solve::theta_value(int num)
{
	switch (num)
	{
	case 0:
	{
		return 0;
	}
	case 1:
	{
		return 0;
	}
	default:
		break;
	}
	return -1;
}

double Solve::S1_gen(double x, double y, double t, int bound1)
{
	switch (bound1)
	{
	case 0:
	{
		//return pow(x, 4);
		return 0;
		//return pow(t, 2);
		break;
	}
	case 1:
	{
		return pow(x, 2);
		//return t * pow(x, 2);
		//return 81 * pow(t, 3);
		break;
	}
	default:
		break;
	}
	return -1;
}

double Solve::distance(Point a, Point b)
{
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

//Генерация правой части
void Solve::gen_vector(vector<double> B, vector<double> B_last, double h, int iter)
{
	vector<double> M_vec_x1(node_count), M_vec_x2(node_count), G_vec(node_count);
	fill(F.begin(), F.end(), 0);

	mult_vec_matrix(M, Q[iter - 1], M_vec_x1);
	mult_vec_matrix(M, Q[iter - 2], M_vec_x2);
	mult_vec_matrix(G, Q[iter - 2], G_vec);

	for (int i = 0; i < node_count; i++)
	{
		F[i] += B[i] / 2 + B_last[i] / 2 + 2 * Arg[iter].hi * M_vec_x1[i] / pow(h, 2) - 
			Arg[iter].hi * M_vec_x2[i] / pow(h, 2) + Arg[iter].sigma * M_vec_x2[i] / (2 * h) - G_vec[i] / 2;
	}
}

//void Solve::vector_addition(vector<double> Vec1, vector<double> Vec2, vector<double>& res, double coeff1, double coeff2)
//{
//	for (int i = 0; i < N; i++)
//	{
//		res[i] += coeff1 * Vec1[i] + coeff2 * Vec2[i];
//	}
//}

//Генерация глобальной матрицы
void Solve::gen_matrix(double h, int iter)
{
	fill(Global.di.begin(), Global.di.end(), 0);
	fill(Global.ggl.begin(), Global.ggl.end(), 0);
	fill(Global.ggu.begin(), Global.ggu.end(), 0);

	for (int i = 0; i < node_count; i++)
	{
		Global.di[i] = G.di[i] / 2 + Arg[iter].hi * M.di[i] / pow(h, 2) + Arg[iter].sigma * M.di[i] / (2 * h);

		for (int j = Global.ig[i]; j < Global.ig[i + 1]; j++)
		{
			Global.ggl[j] = G.ggl[j] / 2 + Arg[iter].hi * M.ggl[j] / pow(h, 2) + Arg[iter].sigma * M.ggl[j] / (2 * h);
			Global.ggu[j] = G.ggu[j] / 2 + Arg[iter].hi * M.ggu[j] / pow(h, 2) + Arg[iter].sigma * M.ggu[j] / (2 * h);
		}
	}
}

void Solve::method()
{
	Matrix Matrix(node_count, triangle_count, G_temp, M_temp);
	vector<double> B_tmp(node_count), B(node_count), B_last(node_count), Q_tmp(node_count, 1);
	double t_now;

	allocate_memory();

	G = Matrix.get_matrix_G();
	M = Matrix.get_matrix_M();
	Unit = Matrix.get_units();

	Global.di.resize(node_count);
	Global.ggl.resize(G.ggl.size());
	Global.ggu.resize(G.ggu.size());
	Global.ig.assign(G.ig.begin(), G.ig.end());
	Global.jg.assign(G.jg.begin(), G.jg.end());

	for (int i = 0; i < node_count; i++)
	{
		Q[0][i] = U_true(Unit[i].x, Unit[i].y, t_start);
		Q[1][i] = U_true(Unit[i].x, Unit[i].y, t_start + t_h);
	}

	t_now = t_start + t_h * 2;

	Matrix.gen_global_vector(t_start);
	B_tmp = Matrix.get_global_vector();

	for (int i = 2; i < time_count; i++)
	{		
		Matrix.gen_global_vector(t_now);
		B = Matrix.get_global_vector();

		if (i % 2 == 0)
		{	
			B_last.assign(B_tmp.begin(), B_tmp.end());
			B_tmp.assign(B.begin(), B.end());
		}
		else
		{
			Matrix.gen_global_vector(t_now - 2 * t_h);
			B_last = Matrix.get_global_vector();
		}

		gen_matrix(t_h, i);
		gen_vector(B, B_last, t_h, i);
		bound2();
		bound1(t_now);
		LOS(Q_tmp);
		Q[i].assign(Q_tmp.begin(), Q_tmp.end());
		t_now += t_h;
	}

}

//Выделение памяти
void Solve::allocate_memory()
{
	Q.resize(time_count);
	U.resize(time_count);
	F.resize(node_count);
	Arg.resize(time_count);
	Unit.resize(node_count);

	for (int i = 0; i < time_count; i++)
	{
		Q[i].resize(node_count);
		U[i].resize(node_count);
	}
}

//LU разложение матрицы
void Solve::LUsq(vector<double>& ggl, vector<double>& ggu, vector<double>& di)
{
	for (int i = 0; i < node_count; i++)
	{
		double sumdi = 0;
		int i0 = Global.ig[i];
		int i1 = Global.ig[i + 1];
		for (int k = i0; k < i1; k++) // к показывает сколько элементов мы обработали в i строке 
		{
			int j = Global.jg[k];// номер столбца  к-го элемента i строки 
			double sumal = 0;
			double sumau = 0;

			int j0 = Global.ig[j];
			int j1 = Global.ig[j + 1];

			int ki = i0;
			int kj = j0;
			for (; ki < k && kj < j1;) // пока есть элементы которые предшествуют к-му элементу
			{
				int j_kl = Global.jg[ki];
				int j_ku = Global.jg[kj];
				if (j_kl == j_ku) // чтобы рассматриваемые элементы не были нулевыми 
				{
					sumal += ggu[kj] * ggl[ki];
					sumau += ggu[ki] * ggl[kj];
					ki++; kj++;
				}
				else if (j_kl < j_ku) ki++;
				else kj++;
			}

			ggu[k] = (ggu[k] - sumau) / di[j];
			ggl[k] = (ggl[k] - sumal) / di[j];
			sumdi += ggl[k] * ggu[k];
		}
		di[i] = sqrt(di[i] - sumdi);
	}
}

//Скалярное произведение векторов
double Solve::scalar_product(vector<double> A1, vector<double> A2)
{
	double sum = 0;

	for (size_t i = 0; i < A1.size(); i++)
	{
		sum += A1[i] * A2[i];
	}

	return sum;
}

//Прямой ход
void Solve::forward(vector<double> ggl, vector<double> di, vector<double> F, vector<double>& X)
{
	fill(X.begin(), X.end(), 0);

	for (int i = 0; i < node_count; i++)
	{
		double sum = 0;

		for (int j = Global.ig[i]; j < Global.ig[i + 1]; j++)
		{
			int st_j = Global.jg[j];
			sum += ggl[j] * X[st_j];
		}
		X[i] = (F[i] - sum) / di[i];
	}
}

//Прямой ход
void Solve::reverse(vector<double> ggu, vector<double> di, vector<double> F, vector<double>& X)
{
	X.assign(F.begin(), F.end());

	for (int i = node_count - 1; i >= 0; i--)
	{
		X[i] /= di[i];

		for (int j = Global.ig[i + 1] - 1; j >= Global.ig[i]; j--)
		{
			int st_j = Global.jg[j];
			X[st_j] -= ggu[j] * X[i];
		}
	}
}

//Умножение матрицы на вектор
void Solve::mult_vec_matrix(GlobalMatrix Mat, vector<double> param, vector<double>& res)
{
	fill(res.begin(), res.end(), 0.0);

	for (int i = 0; i < node_count; i++)
	{
		int i0 = Mat.ig[i];
		int i1 = Mat.ig[i + 1];
		res[i] = Mat.di[i] * param[i];

		for (int j = i0; j < i1; j++)
		{
			int st_j = Mat.jg[j];
			res[i] += param[st_j] * Mat.ggl[j];
			res[st_j] += param[i] * Mat.ggu[j];
		}
	}
}

//Локально оптимальная схема
void Solve::LOS(vector<double>& Q)
{
	double eps = 1e-50, norm, alfa = 0, beta = 0;
	int k = 0, max_iter = 1000;
	vector<double> ggl(Global.ig[node_count]), ggu(Global.ig[node_count]), di(node_count);
	vector<double> mult1(node_count, 0), mult2(node_count), R(node_count), Z(node_count), P(node_count);

	fill(Q.begin(), Q.end(), 1);

	ggl.assign(Global.ggl.begin(), Global.ggl.end());
	ggu.assign(Global.ggu.begin(), Global.ggu.end());
	di.assign(Global.di.begin(), Global.di.end());

	LUsq(ggl, ggu, di);
	mult_vec_matrix(Q, mult1);
	for (int i = 0; i < node_count; i++) mult2[i] = F[i] - mult1[i];

	forward(ggl, di, mult2, R);
	reverse(ggu, di, R, Z);

	mult_vec_matrix(Z, mult1);
	forward(ggl, di, mult1, P);
	norm = sqrt(scalar_product(R, R));

	for (; k < max_iter && norm >= eps; k++)
	{
		double scal_a = scalar_product(P, R);
		double scal_b = scalar_product(P, P);
		alfa = scal_a / scal_b;

		for (int i = 0; i < node_count; i++)
		{
			Q[i] += alfa * Z[i];
			R[i] -= alfa * P[i];
		}
		reverse(ggu, di, R, mult1);
		mult_vec_matrix(mult1, mult2);
		forward(ggl, di, mult2, mult1);

		scal_a = scalar_product(P, mult1);
		beta = -scal_a / scal_b;

		reverse(ggu, di, R, mult2);

		for (int i = 0; i < node_count; i++)
		{
			Z[i] = beta * Z[i] + mult2[i];
			P[i] = mult1[i] + beta * P[i];
		}
		norm = sqrt(scalar_product(R, R));
	}
}

void Solve::mult_vec_matrix(vector<double> param, vector<double>& res)
{
	fill(res.begin(), res.end(), 0.0);

	for (int i = 0; i < node_count; i++)
	{
		int i0 = Global.ig[i];
		int i1 = Global.ig[i + 1];
		res[i] = Global.di[i] * param[i];

		for (int j = i0; j < i1; j++)
		{
			int st_j = Global.jg[j];
			res[i] += param[st_j] * Global.ggl[j];
			res[st_j] += param[i] * Global.ggu[j];
		}
	}
}

//Первое краевое условие
void Solve::bound1(double t_now)
{
	for (auto& edge : S1)
	{
		for (int index : edge.v)
		{
			Global.di[index] = 1.0e+50;
			F[index] = 1.0e+50 * S1_gen(Unit[index].x, Unit[index].y, t_now, edge.valueNo);
		}
	}
}

//Второе краевое условие
void Solve::bound2()
{
	vector<double> tmp_value(4, 0);

	for (auto& edge : S2)
	{
		double h = distance(Unit[edge.v[0]], Unit[edge.v[3]]);

		for (int i = 0; i < 4; i++)
		{
			tmp_value[i] = h * Edge_basis[i] * theta_value(edge.valueNo);
		}

		for (int i = 0; i < 4; i++) F[edge.v[i]] += tmp_value[i];
	}
}

//Чтение данных
void Solve::read_data()
{
	ifstream b1("S1.txt");

	b1 >> boun_num1;
	S1.resize(boun_num1);

	for (int i = 0; i < boun_num1; i++)
	{
		Edge tmp;
		int num;
		for (int j = 0; j < 4; j++)
		{
			b1 >> num;
			tmp.v.push_back(num);
		}
		b1 >> tmp.valueNo;
		S1[i] = tmp;
	}

	b1.close();


	ifstream b2("S2.txt");

	b2 >> boun_num2;
	S2.resize(boun_num2);

	for (int i = 0; i < boun_num2; i++)
	{
		Edge tmp;
		int num;
		for (int j = 0; j < 4; j++)
		{
			b2 >> num;
			tmp.v.push_back(num);
		}
		b2 >> tmp.valueNo;
		S2[i] = tmp;
	}

	b2.close();

	ifstream arguments("arguments.txt");

	double hi_tmp, sigma_tmp;

	arguments >> hi_tmp >> sigma_tmp;

	Arg.resize(time_count);

	for (int i = 0; i < time_count; i++)
	{
		Arg[i].hi = hi_tmp;
		Arg[i].sigma = sigma_tmp;
	}

	arguments.close();
}

vector<vector<double>> Solve::get_U()
{
	double t_now = t_start;

	for (int i = 0; i < time_count; i++)
	{
		for (int j = 0; j < node_count; j++)
		{
			U[i][j] = U_true(Unit[j].x, Unit[j].y, t_now);
		}

		t_now += t_h;
	}

	return U;
}