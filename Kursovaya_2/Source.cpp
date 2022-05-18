#include "Matrix.h"
#include "Solve.h"
#include <iomanip>
#include <vector>

int main()
{
	setlocale(LC_ALL, "Russian");
	int node_count, triangle_count, time_count;
	double t_start, t_end, t_h;
	

	ifstream info("info.txt");
	info >> node_count >> triangle_count >> time_count;
	info.close();

	ifstream time("time.txt");
	time >> t_start >> t_end >> t_h;


	vector<vector<double>> Q, U;

	Solve Met(node_count, triangle_count, time_count, t_start, t_end, t_h);

	Met.method();

	Q = Met.get_Q();
	U = Met.get_U();

	

	double sum = 0;
	double t_now = t_start;

	ofstream out("output.txt");
	out.precision(10);

	for (int i = 0; i < time_count; i++)
	{
		sum = 0;
		out << "time:" << t_now << setw(20) << "Численное значение" << setw(30) << "Точное значение" << setw(20) << "Погрешность" << endl;
		for (int j = 0; j < node_count; j++)
		{
			out << sum << setw(20) << Q[i][j] << setw(30) << U[i][j] << setw(30) << abs(Q[i][j] - U[i][j]) << endl;
			sum++;
		}

		t_now += t_h;
	}
	out.close();
	return 0;
}