/******************************************************************************
 * @details This is a file containing functions regarding quadratures.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/07
 ******************************************************************************/

#include <cassert>
#include <cmath>
#include <functional>

#include <iostream>

namespace quadrature
{
	/******************************************************************************
	 * legendrePolynomial
	 * 
	 * @details    Returns the nth Legendre polynomial.
	 *
	 * @param[in] n 			Gives the degree of the polynomial.
	 * @param[in] x 			The point at which to evaluate the polynomial.
	 * @return 					The Legendre polynomial.
	 ******************************************************************************/
	std::function<double(double)> legendrePolynomial(const int n)
	{
		assert(n >= 0);

		switch(n)
		{
			case 0: return[](double x) -> double
					{
						return 1;
					};
					break;
			case 1: return[](double x) -> double
					{
						return x;
					};
					break;
			default: return[n](double x) -> double
					{
						return common::addFunction  (
														common::constantMultiplyFunction(
																							double(2*n - 1)/n * x,
																							legendrePolynomial(n-1)
																						),
														common::constantMultiplyFunction(
																							-1,
																							common::constantMultiplyFunction(
																																double(n - 1)/n,
																																legendrePolynomial(n-2)
																															)
																						)
													)(x);
					};
		}
	}

	/******************************************************************************
	 * legendrePolynomialDerivative
	 * 
	 * @details    Returns the nth Legendre polynomial's derivative.
	 *
	 * @param[in] n 			Gives the degree of the polynomial.
	 * @param[in] x 			The point at which to evaluate the polynomial.
	 * @return 					The Legendre polynomial derivative.
	 ******************************************************************************/
	std::function<double(double)> legendrePolynomialDerivative(const int n)
	{
		assert(n >= 0);

		switch(n)
		{
			case 0: return[](double x) -> double
					{
						return 0;	
					};
					break;
			case 1: return[](double x) -> double
					{
						return 1;
					};
					break;
			default: return[n](double x) -> double
					{
						return common::addFunction  (
														common::constantMultiplyFunction(
																							double(2*n - 1)/(n - 1) * x,
																							legendrePolynomialDerivative(n-1)
																						),
														common::constantMultiplyFunction(
																							-1,
																							common::constantMultiplyFunction(
																																double(n)/(n - 1),
																																legendrePolynomialDerivative(n-2)
																															)
																						)
													)(x);						
					};

		}
	}

	/******************************************************************************
	 * legendrePolynomialRoots
	 * 
	 * @details    Calculates the ith root of the nth degree Legendre polynomial.
	 *
	 * @param[in] n 			Gives the degree of the polynomial.
	 * @param[in] i 			Which root to return.
	 * @return 					The root of the Legendre polynomial.
	 ******************************************************************************/
	double legendrePolynomialRoot(const int n, const int i)
	{
		double root = -cos(double(2*i + 1)/(2*n)*M_PI);

		while (fabs(legendrePolynomial(n)(root)) >= 1e-5)
			root = root - legendrePolynomial(n)(root)/legendrePolynomialDerivative(n)(root);

		return root;
	}

	/******************************************************************************
	 * legendrePolynomialRoots
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
			roots[i] = legendrePolynomialRoot(n, i);
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
	 * @return 					The quadrature.
	 ******************************************************************************/
	double gaussLegendreQuadrature(const std::function<double(double)> f, const int n)
	{
		double weight, x;
		double h = double(2)/n;
		double quadrature = 0;

		for (int i=0; i<n; ++i)
		{
			x = legendrePolynomialRoot(n, i);
			weight = double(2) / ((1 - pow(x, 2))*pow(legendrePolynomialDerivative(n)(x), 2));
			quadrature += weight * f(x);
		}

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
	 * @return 					The quadrature.
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
}