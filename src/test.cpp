<<<<<<< HEAD
#include "linearSystems.cpp"
#include "quadrature.cpp"
=======
#include "common.cpp"
#include "linearSystems.cpp"
#include "quadrature.cpp"
#include "interpolation.cpp"
>>>>>>> Github/master
#include <cmath>
#include <functional>
#include <iostream>

using namespace std;

<<<<<<< HEAD
=======
double basis0(const double x)
{
	if (0 <= x && x <= 1)
		return 1-x;
	else
		return 0;
}

double basis1(const double x)
{
	if (0 <= x && x <= 1)
		return x;
	else
		return 0;
}

double linearInterpolate(const int n, const double solutions[], const double x)
{
	/*double h = double(2)/n;

	int i = floor(x/h);

	return constantMultiplyFunction(solutions[i], basis1);*/

	assert(-1 <= x && x <= 1);

	double h = double(2)/(n-1);

	bool foundRange = false;
	int i = 0;
	double node1, node2;
	while(!foundRange && i<n-1)
	{
		node1 = -1 + i*h;
		node2 = -1 + (i+1)*h;

		if (node1 <= x && x <= node2)
			foundRange = true;

		++i;
	}

	//function<double(double)> f1 = constantMultiplyFunction(solutions[i], basis0((x - node1)/h));
	//function<double(double)> f2 = constantMultiplyFunction(solutions[i+1], basis1((x - node1)/h));
	function<double(double)> f1 = constantMultiplyFunction(solutions[i], basis0);
	function<double(double)> f2 = constantMultiplyFunction(solutions[i+1], basis1);

	return addFunction(f1, f2)(x);



	
}

>>>>>>> Github/master
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

<<<<<<< HEAD
=======
		/*
>>>>>>> Github/master
		double node;
		double* errors2 = new double[n];
		for (int i=0; i<n; ++i)
		{
			node = -1 + h*i;
<<<<<<< HEAD
			errors2[i] = pow(sol[i] - u(node), 2);
		}

		double integral = trapeziumRule(n, sol, h);
=======
			errors2[i] = pow(sol[i] - u(node), 2); // naughty
		}

		double integral = sqrt(trapeziumRule(n, errors2, h));
		//double integral = trapeziumRule(n, sol, h);
>>>>>>> Github/master

		cout << lastNorm/integral << endl;

		lastNorm = integral;

		delete[] errors2;
<<<<<<< HEAD
=======
		*/

		//cout << linearInterpolate(n, sol, 0.5) << endl;

		//cout << basisL2Norm(n, u);








>>>>>>> Github/master

		delete[] sol;
		delete[] d;
		delete[] c;
		delete[] b;
		delete[] a;
	}

	return 0;
}