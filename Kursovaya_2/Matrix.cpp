#include "Matrix.h"
#include <iomanip>

double Matrix::fun(double x, double y, double t)
{
	//return -2 + 5 * pow(x,2);
	//return -42 * x + 10 * pow(x, 3);
	//return t + 1;
	return -2;
	//return -24 * pow(x, 2) + pow(x, 4);
	//return 3 * pow(x, 2) - 2*t;
	//return 3 * pow(x, 4) * pow(t, 2) + 6 * pow(x, 4) * t - 12 * pow(x, 2) * pow(t, 3);
	//return 2 * t * pow(x, 3) + 2 * pow(x, 3) - 6 * pow(t, 2) * x;
}

//Чтение данных
void Matrix::read_matrix()
{
	ifstream xy("xy.txt");

	for (int i = 0; i < node_count; i++)
	{
		double x, y;
		xy >> x >> y;
		Point Temp(x, y);
		Unit.push_back(Temp);
	}
	xy.close();

	ifstream arg("lambda.txt");

	double lambda;
	arg >> lambda;

	for (int i = 0; i < triangle_count; i++)
	{
		arg_lambda[i] = lambda;
	}

	arg.close();

	ifstream tri("triangle.txt");

	for (int i = 0; i < triangle_count; i++)
	{
		vector<int> temp(10), temp_sort(10), new_sort_temp(10);
		for (int j = 0; j < 10; j++)
		{
			tri >> temp[j];
			temp_sort[j] = temp[j];
		}
		sort(temp_sort.begin(), temp_sort.end());

		for (int j = 0; j < 10; j++)
		{
			for (int k = 0; k < 10; k++)
			{
				if (temp_sort[j] == temp[k])
				{
					new_sort_temp[j] = k;
				}
			}
		}

		Triangle T(temp, new_sort_temp);
		Ind_tri.push_back(T);
	}
	tri.close();
}

//Выделение памяти
void Matrix::memory_allocation()
{
	Unit.reserve(node_count);
	Ind_tri.reserve(triangle_count);
	arg_lambda.resize(triangle_count);
	M_glob.ig.resize(node_count + 1);
	M_glob.di.resize(node_count);
	G_glob.ig.resize(node_count + 1);
	G_glob.di.resize(node_count);
	M_coeff.resize(10);
	G_coeff.resize(10);
	G_local.resize(10);
	F_global.resize(node_count);
	F_local.resize(10);
	Alfa_vec.resize(3);

	for (int i = 0; i < 10; i++)
	{
		M_coeff[i].resize(i + 1);
		G_coeff[i].resize(i + 1);
		G_local[i].resize(10);
	}
}

//Генерация гломальной матрицы G
void Matrix::gen_global_G()
{
	gen_G();

	for (int i = 0; i < triangle_count; i++)
	{
		gen_local_G(Ind_tri[i], arg_lambda[i]);
		add_to_glob(Ind_tri[i], G_glob, G_local);
	}
}

//Генерация глобальной матрицы M
void Matrix::gen_global_M()
{
	M_glob.ig.assign(G_glob.ig.begin(), G_glob.ig.end());
	M_glob.jg.assign(G_glob.jg.begin(), G_glob.jg.end());
	gen_M();

	for (int i = 0; i < triangle_count; i++)
	{
		add_to_glob(Ind_tri[i], M_glob, M_coeff);
	}
}

//Генерация локальной матрицы G
void Matrix::gen_local_G(Triangle& tri, double lambda)
{
	vector<Point> Temp(3), LocalTri(10);
	double D, sumG;

	for (int i = 0; i < 10; i++)
	{
		fill(G_local[i].begin(), G_local[i].end(), 0);
	}

	for (int i = 0; i < 10; i++)
	{
		LocalTri[i] = Unit[tri.Index_first[i]];
	}

	for (int i = 0; i < 3; i++) Temp[i] = LocalTri[i];

	D = detD(Temp);
	gen_Alfa(Temp, D);

	for (size_t i = 0; i < 10; i++)
	{
		for (size_t j = 0; j <= i; j++)
		{
			sumG = 0;

			for (size_t k = 0; k < G_coeff[i][j].size(); k++)
			{
				int gradno1 = G_coeff[i][j][k].grad1;
				int gradno2 = G_coeff[i][j][k].grad2;
				sumG += (Alfa_vec[gradno1].alfa1 * Alfa_vec[gradno2].alfa1 + Alfa_vec[gradno1].alfa2 * Alfa_vec[gradno2].alfa2) * G_coeff[i][j][k].coeff;
			}

			G_local[i][j] = lambda * sumG;
			G_local[j][i] = G_local[i][j];
		}
	}

}

//Добавление локальной матрицы в глобальную
void Matrix::add_to_glob(Triangle& tri, GlobalMatrix& Glob_mat, vector<vector<double>> Local_mat)
{
	for (int i = 0; i < 10; i++)
	{
		vector<Point> Temp(3);
		double D;

		int diElem = tri.Index_first[i];
		
		for (int j = 0; j < 3; j++)
		{
			Temp[j] = Unit[tri.Index_first[j]];
		}

		D = detD(Temp);

		Glob_mat.di[diElem] += Local_mat[i][i] * abs(D);

		for (int j = 0; j < i; j++)
		{

			int a = tri.Index_first[i];
			int b = tri.Index_first[j];
			if (a < b) swap(a, b);

			auto begin = Glob_mat.jg.begin() + Glob_mat.ig[a];
			if (Glob_mat.ig[a + 1] > Glob_mat.ig[a])
			{
				auto end = Glob_mat.jg.begin() + Glob_mat.ig[a + 1] - 1;
				auto iter = lower_bound(begin, end, b);
				auto index = iter - Glob_mat.jg.begin();
				Glob_mat.ggl[index] += Local_mat[i][j] * abs(D);
				Glob_mat.ggu[index] += Local_mat[i][j] * abs(D);
			}
		}
	}
}

//Построение портрета глобальной матрицы
void Matrix::portrait()
{
	vector<set<int>> Beg(node_count);

	for (auto& ielem : Ind_tri)
	{
		for (int i = 1; i < 10; i++)
		{
			for (int j = 0; j < i; j++)
			{
				auto a = ielem.Index_first[i];
				auto b = ielem.Index_first[j];
				if (a < b) swap(a, b);

				Beg[a].insert(b);
			}
		}
	}

	G_glob.ig[0] = G_glob.ig[1] = 0;

	for (int i = 2; i <= node_count; i++)
	{
		int col = G_glob.ig[i - 1];
		G_glob.ig[i] = col + Beg[i - 1].size();
	}

	G_glob.jg.resize(G_glob.ig[node_count]);
	G_glob.ggl.resize(G_glob.ig[node_count]);
	G_glob.ggu.resize(G_glob.ig[node_count]);
	M_glob.jg.resize(G_glob.ig[node_count]);
	M_glob.ggl.resize(G_glob.ig[node_count]);
	M_glob.ggu.resize(G_glob.ig[node_count]);

	for (int i = 1, k = 0; i < node_count; i++)
	{
		for (int j : Beg[i])
		{
			G_glob.jg[k] = j;
			k++;
		}
	}
}

//Генерация коэффициентов альфа для локаьной матрицы
void Matrix::gen_Alfa(vector<Point>& P, double det)
{
	Alfa_vec[0].alfa1 = (P[1].x - P[2].x) / det;
	Alfa_vec[0].alfa2 = (P[2].y - P[1].y) / det;
	Alfa_vec[1].alfa1 = (P[2].x - P[0].x) / det;
	Alfa_vec[1].alfa2 = (P[0].y - P[2].y) / det;
	Alfa_vec[2].alfa1 = (P[0].x - P[1].x) / det;
	Alfa_vec[2].alfa2 = (P[1].y - P[0].y) / det;
}

//Определитель матрицы D
double Matrix::detD(vector<Point>& P)
{
	double D;
	D = (P[1].x - P[0].x) * (P[2].y - P[0].y) - (P[2].x - P[0].x) * (P[1].y - P[0].y);
	return D;
}

//Генерация матрицы массы
void Matrix::gen_M()
{
	double coeff = 0;

	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			coeff = 0;
			for (auto& x : M_temp[i])
			{
				for (auto& y : M_temp[j])
				{
					int v1 = x.v1 + y.v1;
					int v2 = x.v2 + y.v2;
					int v3 = x.v3 + y.v3;

					coeff += x.coeff * y.coeff * factorial(v1) * factorial(v2) * factorial(v3) /
						factorial(v1 + v2 + v3 + 2);

				}
			}
			M_coeff[i][j] = coeff;
		}
	}
}

//Гереация матрицы жесткости
void Matrix::gen_G()
{
	double coeff;

	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			for (auto& x : G_temp[i])
			{
				for (auto& y : G_temp[j])
				{
					int v1 = x.v1 + y.v1;
					int v2 = x.v2 + y.v2;
					int v3 = x.v3 + y.v3;

					coeff = x.coeff * y.coeff * factorial(v1) * factorial(v2) * factorial(v3) /
						factorial(v1 + v2 + v3 + 2);
					G_coeff[i][j].push_back(GradMult(x.gradNo, y.gradNo, coeff));
				}
			}
		}
	}
}

void Matrix::gen_global_vector(double t_now)
{

	fill(F_global.begin(), F_global.end(), 0);

	for (int i = 0; i < triangle_count; i++)
	{
		gen_local_vector(Ind_tri[i], t_now);

		for (int j = 0; j < 10; j++)
		{
			int diElem = Ind_tri[i].Index_first[j];

			F_global[diElem] += F_local[j];
		}
	}
}

//Генерация локального вектора
void Matrix::gen_local_vector(Triangle& tri, double t_now)
{
	vector<double> p1 =
	{ 0.0451890097844, 0.0451890097844, 0.9096219804312, 0.7475124727339, 0.2220631655373, 0.7475124727339,
		0.2220631655373, 0.0304243617288, 0.0304243617288, 0.1369912012649, 0.6447187277637, 0.1369912012649, 0.2182900709714,
		0.2182900709714, 0.6447187277637, 0.0369603304334, 0.4815198347833, 0.4815198347833, 0.4036039798179, 0.4036039798179,
		0.1927920403641
	};
	vector<double> p2 =
	{ 0.0451890097844, 0.9096219804312, 0.0451890097844, 0.0304243617288, 0.0304243617288, 0.2220631655373,
		0.7475124727339, 0.7475124727339, 0.2220631655373, 0.2182900709714, 0.2182900709714, 0.6447187277637, 0.6447187277637,
		0.1369912012649, 0.1369912012649, 0.4815198347833, 0.0369603304334, 0.4815198347833, 0.1927920403641, 0.4036039798179,
		0.4036039798179
	};
	// веса
	vector<double> w =
	{ 0.0519871420646, 0.0519871420646, 0.0519871420646, 0.0707034101784, 0.0707034101784, 0.0707034101784,
		0.0707034101784, 0.0707034101784, 0.0707034101784, 0.0909390760952, 0.0909390760952, 0.0909390760952, 0.0909390760952,
		0.0909390760952, 0.0909390760952, 0.1032344051380, 0.1032344051380, 0.1032344051380, 0.1881601469167, 0.1881601469167,
		0.1881601469167
	};

	// локальные координаты для каждой квадратуры
	vector<double> LQ(3);
	vector<Point> node(3);
	fill(F_local.begin(), F_local.end(), 0);
	for (int i = 0; i < 3; i++) node[i] = Unit[tri.Index_first[i]];
	// значения x и y в каждой квадратуре
	double xQ = 0, yQ = 0;
	// вспомогательная переменная
	double s = 0;

	double D = abs(detD(node));

	for (int i = 0; i < 21; i++)
	{
		LQ[0] = p1[i];
		LQ[1] = p2[i];
		LQ[2] = 1 - LQ[1] - LQ[0];

		xQ = node[0].x * LQ[0] + node[1].x * LQ[1] + node[2].x * LQ[2];
		yQ = node[0].y * LQ[0] + node[1].y * LQ[1] + node[2].y * LQ[2];

		s = 0.25 * w[i] * fun(xQ, yQ, t_now);

		vector<double> valuePsi(10, 0);

		for (int j = 0; j < 10; j++)
		{
			for (auto& elem : M_temp[j])
			{
				valuePsi[j] += elem.coeff * pow(LQ[0], elem.v1) * pow(LQ[1], elem.v2) * pow(LQ[2], elem.v3);
			}
		}


		for (int j = 0; j < 10; j++)
			F_local[j] += valuePsi[j] * s * D;
	}
}

double Matrix::factorial(int v)
{
	double res = 1;
	for (int i = 2; i <= v; i++)
	{
		res *= i;
	}
	return res;
}
