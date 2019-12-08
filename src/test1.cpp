#include "common.cpp"
#include "linearSystems.cpp"
#include "quadrature.cpp"
#include "interpolation.cpp"
#include <cmath>
#include <functional>
#include <iostream>

using namespace std;
using namespace common;

double basis1(const double x)
{
	if (-1 <= x && x <= 1)
		return (1-x)/2;
	else
		return 0;
}

double basis2(const double x)
{
	if (-1 <= x && x <= 1)
		return (1+x)/2;
	else
		return 0;
}

double linearInterpolate(const int n, const double solutions[], const double x)
{
	assert(-1 <= x && x <= 1);

	double h = double(2)/(n-1);

	bool foundRange = false;
	int i = 0;
	double node1, node2;
	do 
	{
		node1 = -1 + i*h;
		node2 = -1 + (i+1)*h;

		if (node1 <= x && x <= node2)
			foundRange = true;
		else
			++i;
	} while(!foundRange && i<n);

	function<double(double)> f1 = constantMultiplyFunction(solutions[i], basis1);
	function<double(double)> f2 = constantMultiplyFunction(solutions[i+1], basis2);

	return addFunction(f1, f2)(2*(x - node1)/(node2 - node1) - 1);
}

double u(double x)
{
	return (x-1)*(x+1)/2;
}

double solutionL2Norm(const int n, const double solution[], function<double(double)> u)
{
	double norm = 0;

	double node1, node2;
	double h = double(2)/(n-1);
	function<double(double)> u_;
	function<double(double)> error;
	function<double(double)> error2;

	// Loops over each element.
	for (int i=0; i<n-1; ++i)
	{
		node1 = -1 + h*i;
		node2 = node1 + h;
		u_ = addFunction(constantMultiplyFunction(solution[i], basis1), constantMultiplyFunction(solution[i+1], basis2));
		error = addFunction(transformFunction(u, node1, node2), constantMultiplyFunction(-1, u_));
		error2 = multiplyFunction(error, error);
		norm += gaussLegendreQuadrature(error2, 4)*(node2-node1)/2;
		//cout << "x = " << -1+h*i << " " << error(-1) << " " << error(0) << " " << error(1) << endl;
		//for (int j=0; j<100; ++j)
			//cout << "WOW" << j << " " << solution[i] << " " << solution[i+1] << " " << u_(-1+j*double(2)/99) << " " << u(-1+h*i+j*h*double(1)/(99)) << endl;
			//cout << "WOW" << j << " " << solution[i] << " " << solution[i+1] << " " << u_(-1+j*double(2)/99) << " " << transformFunction(u, -1+h*i, -1+h*(i+1))(-1+j*double(2)/99) << endl;
		//cout << "WOW: " << gaussLegendreQuadrature(error2, n)*(node2-node1)/2 << endl;
	}

	return sqrt(norm);
}

double f(double x)
{
	return 1;
}

int main()
{
	double lastNorm = 0;

	//for (int n=2; n<=32; n*=2)
	for (int n=2; n<=1024; n*=2)
	//for (int n=2; n<=8; n*=2)
	//for (int n=4; n==4; n*=2)
	{
		double h = double(2)/(n+1);

		//function<double(double)> f = F;

		double* a   = new double[n-1];
		double* b   = new double[n];
		double* c   = new double[n-1];
		double* d   = new double[n];
		double* sol = new double[n+2];

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

		sol[0] = 0;
		sol[n+1] = 0;

		thomasInvert(n, a, b, c, d, sol+1);

		//for (int i=0; i<n+2; ++i)
			//cout << sol[i] << endl;

		cout << endl;

		double norm = solutionL2Norm(n+2, sol, u);

		cout << norm << endl;

		cout << lastNorm/norm << endl;

		lastNorm = norm;

		delete[] sol;
		delete[] d;
		delete[] c;
		delete[] b;
		delete[] a;
	}

	cout << "All done!" << endl;

	return 0;
}