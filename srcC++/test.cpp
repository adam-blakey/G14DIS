#include "linearSystems.cpp"
#include "quadrature.cpp"
#include <cmath>
#include <functional>
#include <iostream>

using namespace std;

double u(double x)
{
	return x*(1-pow(x, 4));
}

double F(double x)
{
	return -20*pow(x, 3);
}

int main()
{
	double lastNorm = 0;

	for (int i=0; i<12; ++i)
	{
		int n = pow(2, i+2);
		double h = double(2)/(n-1);

		function<double(double)> f = F;

		double* a   = new double[n-1];
		double* b   = new double[n];
		double* c   = new double[n-1];
		double* d   = new double[n];
		double* sol = new double[n];

		for (int i=0; i<n-1; ++i)
		{
			a[i] = 1;
			c[i] = 1;
		}
		for (int i=0; i<n; ++i)
		{
			b[i] = -2;
		}
		for (int i=0; i<n; ++i)
		{
			d[i] = pow(h, 2)*f(-1 + (i+1)*h);
		}

		thomasInvert(n, a, b, c, d, sol);

		double node;
		double* errors2 = new double[n];
		for (int i=0; i<n; ++i)
		{
			node = -1 + h*i;
			errors2[i] = pow(sol[i] - u(node), 2);
		}

		double integral = trapeziumRule(n, sol, h);

		cout << lastNorm/integral << endl;

		lastNorm = integral;

		delete[] errors2;

		delete[] sol;
		delete[] d;
		delete[] c;
		delete[] b;
		delete[] a;
	}

	return 0;
}