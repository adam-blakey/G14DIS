/*=============================================================================
 * QUADRATURE.CPP
 *
 * This is a file containing functions regarding quadratures.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/11/4
  =============================================================================*/

#include <cassert>
#include <cmath>
#include <functional>

#include <iostream>

/******************************************************************************
 * legendrePolynomial
 * 
 * @details    Returns the nth Legendre polynomial.
 *
 * @param[in] n 			Gives the degree of the polynomial.
 * @param[in] x 			The point at which to evaluate the polynomial.
 ******************************************************************************/

double legendrePolynomial(const double x, const int n)
{
	assert(n >= 0);

	if (n==0)
		return 1;
	else if (n==1)
		return x;
	else
		return double(2*n - 1)/n * x * legendrePolynomial(x, n-1) - double(n - 1)/n * legendrePolynomial(x, n-2);
}

/******************************************************************************
 * legendrePolynomialDerivative
 * 
 * @details    Returns the nth Legendre polynomial's derivative.
 *
 * @param[in] n 			Gives the degree of the polynomial.
 * @param[in] x 			The point at which to evaluate the polynomial.
 ******************************************************************************/

double legendrePolynomialDerivative(const double x, const int n)
{
	assert(n >= 0);

	if (n==0)
		return 0;
	else if (n==1)
		return 1;
	else
		return double(2*n - 1)/(n - 1) * x * legendrePolynomialDerivative(x, n-1)
				- double(n)/(n - 1) * legendrePolynomialDerivative(x, n-2);
}

/******************************************************************************
 * legendrePolynomialRoot
 * 
 * @details    Calculates the roots of the real roots of the nth Legendre
 * 				polynomial at the point x.
 *
 * @param[in] n 			Gives the degree of the polynomial.
 * @param[out] roots 		The n roots of the polynomial.
 ******************************************************************************/

void legendrePolynomialRoots(const int n, double roots[])
{
	for (int i=0; i<n; ++i)
	{
		roots[i] = -cos(double(2*i - 1)/(2*n)*M_PI);

		while (fabs(legendrePolynomial(roots[i], n)) >= 1e-3)
		{
			roots[i] = roots[i] - legendrePolynomial(roots[i], n)/legendrePolynomialDerivative(roots[i], n);
		}
	}		
}

/******************************************************************************
 * gaussLegendreQuadrature
 * 
 * @details    Calculates the Gauss-Legendre quadrature of f between -1 and 1.
 *
 * @param[in] f 			A function, f, to integrate with.
 * @param[in] n 			What degree of Legendre polynomial to approximate
 * 							 with.
 ******************************************************************************/

double gaussLegendreQuadrature(const std::function<double(double)> f, const int n)
{
	double weight, x;
	double h = double(2)/n;
	double quadrature = 0;

	double* roots = new double[n];

	legendrePolynomialRoots(n, roots);

	for (int i=0; i<n; ++i)
	{
		x = roots[i];
		weight = double(2) / ((1 - pow(x, 2))*pow(legendrePolynomialDerivative(x, n), 2));
		quadrature += weight * f(x);
		//std::cout << "QUAD: " << weight * f(x) << std::endl;
	}

	delete[] roots;

	return quadrature;
}

/******************************************************************************
 * trapeziumRule
 * 
 * @details    Uses composite trapezium rule with n-1 equal strips to give a
 * 				quadrature.
 *
 * @param[in] n 			How many values we're approximating the integral
 * 							 with.
 * @param[in] fValues 		The function values at the sides of each strip.
 * @param[in] h 			The trip width.
 ******************************************************************************/

double trapeziumRule(const int n, const double fValues[], const double h)
{
	double answer = 0;

	answer += fValues[0]/2;
	for (int i=1; i<n-1; ++i)
	{
		answer += fValues[i];
	}
	answer += fValues[n-1]/2;

	return h/(n-1) * answer;
}