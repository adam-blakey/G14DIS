/*=============================================================================
 * COMMON.CPP
 *
 * This is a file containing functions regarding common functions.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/11/11
  =============================================================================*/

#include <cmath>
#include <functional>

using namespace std;

void setToZero(const int n, double array[])
{
	for (int i=0; i<n; ++i)
	{
		array[i] = 0;
	}
}

template <class F, class G>
decltype(auto) composeFunction(F&& f, G&& g)
{
	return [=](auto x) { return f(g(x)); };
}

function<double(double)> transformFunction(function<double(double)> f, double xj, double xjp1)
{
	//return [=](double x)->double{ return f(((x + 1)/2 + xj)*(xjp1 - xj)); };
	return [=](double x)->double{ return f((x+1)*(xjp1-xj)/2 + xj); };
}

function<double(double)> addFunction(function<double(double)> f, function<double(double)> g)
{
	return [=](double x)->double{ return f(x) + g(x); };
}

/*template <class F, class G>
decltype(auto) multiplyFunction(F&& f, G&& g)
{
	return [&](auto x) { return f(x) * g(x); };
}*/

function<double(double)> multiplyFunction(function<double(double)> f, function<double(double)> g)
{
	return [=](double x)->double{ return f(x) * g(x); };
}

function<double(double)> constantMultiplyFunction(double a, function<double(double)> f)
{
	return [=](double x)->double{ return a * f(x); };
}

double basisL2Norm(const int n, std::function<double(double)> u)
{
	double norm = 0;
	double h = double(2)/(n-1);
	double node;

	for (int i=0; i<n; ++i)
	{
		node = -1 + i*h;
		norm += pow(u(node), 2)*h;
	}

	return sqrt(norm);
}

void tridiagonalVectorMultiplication(const int n, const double a[], const double b[], const double c[], const double d[], double solution[])
{
	for (int i=0; i<n; ++i)
	{
		if (i==0)
			solution[i] = b[0]*d[0] + c[0]*d[1];
		else if (i==n-1)
			solution[i] = a[n-2]*d[n-2] + b[n-1]*d[n-1];
		else
			solution[i] = a[i-1]*d[i-1] + b[i]*d[i] + c[i]*d[i+1];
	}
}