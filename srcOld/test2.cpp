#include "common.cpp"
#include "linearSystems.cpp"
#include "quadrature.cpp"
#include "interpolation.cpp"
#include <cmath>
#include <functional>
#include <iostream>

using namespace std;

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

double solutionL2Norm(const int n, const double solution[], function<double(double)> u)
{
	double norm = 0;

	function<double(double)> u_;
	function<double(double)> error;
	function<double(double)> error2;

	// Loops over each element.
	for (int i=0; i<n-1; ++i)
	{
		u_ = addFunction(constantMultiplyFunction(solution[i], basis1), constantMultiplyFunction(solution[i+1], basis2));
		error = addFunction(u, constantMultiplyFunction(-1, u_));
		error2 = multiplyFunction(error, error);
		norm += gaussLegendreQuadrature(error2, n);
		//cout << "WOW: " << gaussLegendreQuadrature(error2, n) << endl;
	}

	return sqrt(norm);
}

double F(double x)
{
	return 1;
}

double u(double x)
{
	return (x-1)*(x+1);
}

int main()
{
	double lastNorm = 0;
	for (int n=1; n<=16; n*=2)
	{
		double* sols = new double[n];
		for (int i=0; i<n; ++i)
			sols[i] = 1 + double(1)/pow(n, 2);
	
		double norm = solutionL2Norm(n, sols, F);
	
		cout << lastNorm/norm << endl;

		lastNorm = norm;
	}

	return 0;
}