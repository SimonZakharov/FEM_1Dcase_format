#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctime>

using namespace std;

//	�������� �������� (����� � ������ �����)
const double EPS = 1e-10;

//	����� ����� �����
long n = 0;

//	�������, ������� � ������ �����
double F(double x)
{
	return -9.03 * exp(0.1 * x);
}

//	������� ������� - ��� ��������
double U(double x)
{
	return exp(0.1 * x);
}

//	����� ��������
void TMA(double *diag, double *up, double *down, double *f, double *u)
{
	double *d, *l;
	d = new double[n]; l = new double[n];
	//	������ ��������
	d[0] = -up[0] / diag[0];
	l[0] = f[0] / diag[0];
	for (int i = 1; i < n; i++)
	{
		if (i != n - 1)
			d[i] = -up[i] / (diag[i] + down[i - 1] * d[i - 1]);
		else d[i] = 0;
		l[i] = (f[i] - down[i - 1] * l[i - 1]) / (diag[i] + down[i - 1] * d[i - 1]);
	}
	//	�������� ��������
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
	//	���������� ������ � ����� ������� + ����������� ���������� ���� �����
	double p0 = 0, c1 = 0, p1 = 0, coeff = 0.;
	//	���������� ��� �������� ���� �����
	double h, h0;
	//	����� ���������
	long m = 0;
	//	������ ��������� ����� �����
	double *x;
	//	���������� ������� ���������
	double G[2][2];
	//	���������� ������� ����
	double M[2][2];
	//	���������� ������ ������ �����
	double b[2];
	//	������� �������
	double *diag;	//	������������ ��������
	double *up;		//	������������ ��� ����������
	double *down;		//	������������ ��� ����������
	//	���������� ������ ������ �����
	double *f;
	//	��������� � ���������
	double lambda = 0, gamma = 0;
	//	������ �������
	double *u;
	//	����������, ����������� ��� ����� ��������� �������
	double mC = 0, nC = 0;
	//	��������� ������
	ifstream input; input.open("input.txt");
	input >> p0 >> p1 >> h0 >> coeff;
	input >> lambda >> gamma;
	input >> mC >> nC;
	input.close();
	c1 = p0;
	h = h0;
	//	�������� ������������� �����
	while (c1 < p1)
	{
		n++;
		c1 += h;
		h *= coeff;
	}
	n++;
	//	������� ������ ��� ����������
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
	//	������� �����
	ofstream output; output.open("output.txt");

	output << "������� ���������: " << -lambda << " * d^2u / dx^2 + " << gamma << " * u = -21 * e^(2x);\n";
	output << "����� �������������.\n��������� ��� ����� h0 = " << h0 << "\n" << "����������� ���������� k = " << coeff << endl;
	output << "��������� �������: [" << p0 << "; " << p1 << "]\n";
	output << "�� ����� ������� ��������� ������� ������ ��������� ������� ������� (�������), �� ������ ������� � ��������� ������� ������� (������������).\n\n";
	output << "������� ������� ���������� ������� ��������.\n";

	output << "���� �����:\n\n";
	for (int i = 0; i < n; i++)
		output << i << "\t" << x[i] << endl;
	m = n - 1;
	//	���� ������ �� ���������
	for (int i = 0; i < m; i++)
	{
		h = x[i + 1] - x[i];
		//	��������� ���������� ������� ���������
		G[0][0] = G[1][1] = lambda / h;
		G[0][1] = G[1][0] = -lambda / h;
		//	��������� ��������� ������� ����
		M[0][0] = M[1][1] = gamma * h / 3;
		M[0][1] = M[1][0] = gamma * h / 6;
		//	��������� ���������� ������ ������ �����
		b[0] = h * (2 * F(x[i]) + F(x[i + 1])) / 6;
		b[1] = h * (F(x[i]) + 2 * F(x[i + 1])) / 6;
		//	�������� ���������� ������� �������
		diag[i] += G[0][0] + M[0][0];
		diag[i + 1] += G[1][1] + M[1][1];
		up[i] += G[0][1] + M[0][1];
		down[i] += G[1][0] + M[1][0];
		//	�������� ���������� ������ ������ �����
		f[i] += b[0]; f[i + 1] += b[1];
	}
	//	����� ������������ ��������� �������
	//f[n - 1] += nC;
	f[n - 1] = U(x[n - 1]);
	diag[n - 1] = 1; down[n - 2] = 0;
	//	����� ������� ��������� �������
	f[0] = U(x[0]);
	diag[0] = 1; up[0] = 0;
	//	������ ���������� ������ ���������� ������� � ������� ���������
	clock_t time_tma;
	time_tma = clock();
	TMA(diag, up, down, f, u);
	time_tma = clock() - time_tma;
	//	����� �������
	output << "\n������ �������:\n\n";
	for (int i = 0; i < n; i++)
	{
		output << U(x[i]) << "\t" << u[i]  << "\t" << abs(U(x[i]) - u[i]) << endl;
	}
	//	������ ������������� ����������� �������
	double e1 = 0, e2 = 0;
	for (int i = 0; i < n; i++)
	{
		e1 += (U(x[i]) - u[i]) * (U(x[i]) - u[i]);
		e2 += U(x[i]) * U(x[i]);
	}
	output << "����� ������������� �����������:\n";
	output << e1 / e2 << endl;
	output << "����� ������ ��������: " << time_tma << "ms";
	output.close();
	return EXIT_SUCCESS;
}