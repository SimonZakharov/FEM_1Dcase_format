#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctime>

using namespace std;

//	желаемая точность (нужна в методе Якоби)
const double EPS = 1e-10;

//	число узлов сетки
long n = 0;

//	функция, стоящая в правой части
double F(double x)
{
	return -9.03 * exp(0.1 * x);
}

//	функция решения - для проверки
double U(double x)
{
	return exp(0.1 * x);
}

//	метод прогонки
void TMA(double *diag, double *up, double *down, double *f, double *u)
{
	double *d, *l;
	d = new double[n]; l = new double[n];
	//	прямая прогонка
	d[0] = -up[0] / diag[0];
	l[0] = f[0] / diag[0];
	for (int i = 1; i < n; i++)
	{
		if (i != n - 1)
			d[i] = -up[i] / (diag[i] + down[i - 1] * d[i - 1]);
		else d[i] = 0;
		l[i] = (f[i] - down[i - 1] * l[i - 1]) / (diag[i] + down[i - 1] * d[i - 1]);
	}
	//	обратная прогонка
	u[n - 1] = l[n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		u[i] = d[i] * u[i + 1] + l[i];
	}
	delete []d;
	delete []l;
}

int main(int argc, char *argv[])
{
	n = 0;
	//	координаты начала и конца отрезка + коэффициент увеличения шага сетки
	double p0 = 0, c1 = 0, p1 = 0, coeff = 0.;
	//	переменные для хранения шага сетки
	double h, h0;
	//	число элементов
	long m = 0;
	//	массив координат узлов сетки
	double *x;
	//	элементная матрица жесткости
	double G[2][2];
	//	элементная матрица масс
	double M[2][2];
	//	элементный вектор правой части
	double b[2];
	//	матрица системы
	double *diag;	//	диагональные элементы
	double *up;		//	коэффициенты над диагональю
	double *down;		//	коэффициенты под диагональю
	//	глобальный вектор правой части
	double *f;
	//	параметры в уравнении
	double lambda = 0, gamma = 0;
	//	вектор решения
	double *u;
	//	переменные, необходимые для учета граничных условий
	double mC = 0, nC = 0;
	//	считываем данные
	ifstream input; input.open("input.txt");
	input >> p0 >> p1 >> h0 >> coeff;
	input >> lambda >> gamma;
	input >> mC >> nC;
	input.close();
	c1 = p0;
	h = h0;
	//	построим неравномерную сетку
	while (c1 < p1)
	{
		n++;
		c1 += h;
		h *= coeff;
	}
	n++;
	//	выделим память под переменные
	x = new double[n];
	f = new double[n]; for (int i = 0; i < n; i++) f[i] = 0;
	u = new double[n];
	diag = new double[n]; for (int i = 0; i < n; i++) diag[i] = 0;
	up = new double[n]; for (int i = 0; i < n; i++) up[i] = 0;
	down = new double[n]; for (int i = 0; i < n; i++) down[i] = 0;
	c1 = p0; h = h0;
	for (int i = 0; i < n - 1; i++)
	{
		x[i] = c1;
		c1 += h;
		h *= coeff;
	}
	x[n - 1] = p1;
	//	выведем сетку
	ofstream output; output.open("output.txt");

	output << "Решение уравнения: " << -lambda << " * d^2u / dx^2 + " << gamma << " * u = -21 * e^(2x);\n";
	output << "Сетка неравномерная.\nНачальный шаг сетки h0 = " << h0 << "\n" << "Коэффициент разрежения k = " << coeff << endl;
	output << "Расчетная область: [" << p0 << "; " << p1 << "]\n";
	output << "На левой границе расчетной области задано граничное условие Дирихле (главное), на правой границе – граничное условие Неймана (естественное).\n\n";
	output << "Решение системы произведем методом прогонки.\n";

	output << "Узлы сетки:\n\n";
	for (int i = 0; i < n; i++)
		output << i << "\t" << x[i] << endl;
	m = n - 1;
	//	идем циклом по элементам
	for (int i = 0; i < m; i++)
	{
		h = x[i + 1] - x[i];
		//	формируем элементную матрицу жесткости
		G[0][0] = G[1][1] = lambda / h;
		G[0][1] = G[1][0] = -lambda / h;
		//	формируем локальную матрицу масс
		M[0][0] = M[1][1] = gamma * h / 3;
		M[0][1] = M[1][0] = gamma * h / 6;
		//	формируем элементный вектор правой части
		b[0] = h * (2 * F(x[i]) + F(x[i + 1])) / 6;
		b[1] = h * (F(x[i]) + 2 * F(x[i + 1])) / 6;
		//	собираем глобальную матрицу системы
		diag[i] += G[0][0] + M[0][0];
		diag[i + 1] += G[1][1] + M[1][1];
		up[i] += G[0][1] + M[0][1];
		down[i] += G[1][0] + M[1][0];
		//	собираем глобальный вектор правой части
		f[i] += b[0]; f[i + 1] += b[1];
	}
	//	учтем естественное граничное условие
	//f[n - 1] += nC;
	f[n - 1] = U(x[n - 1]);
	diag[n - 1] = 1; down[n - 2] = 0;
	//	учтем главное граничное условие
	f[0] = U(x[0]);
	diag[0] = 1; up[0] = 0;
	//	теперь необходимо решить полученную систему и вывести результат
	clock_t time_tma;
	time_tma = clock();
	TMA(diag, up, down, f, u);
	time_tma = clock() - time_tma;
	//	вывод решения
	output << "\nВектор решения:\n\n";
	for (int i = 0; i < n; i++)
	{
		output << U(x[i]) << "\t" << u[i]  << "\t" << abs(U(x[i]) - u[i]) << endl;
	}
	//	оценим относительную погрешность решения
	double e1 = 0, e2 = 0;
	for (int i = 0; i < n; i++)
	{
		e1 += (U(x[i]) - u[i]) * (U(x[i]) - u[i]);
		e2 += U(x[i]) * U(x[i]);
	}
	output << "Норма относительной погрешности:\n";
	output << e1 / e2 << endl;
	output << "Время работы решателя: " << time_tma << "ms";
	output.close();
	return EXIT_SUCCESS;
}