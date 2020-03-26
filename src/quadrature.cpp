/******************************************************************************
 * @details This is a file containing functions regarding quadratures.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/07
 ******************************************************************************/
#include "common.hpp"
#include "quadrature.hpp"
#include <cassert>
#include <cmath>
#include <functional>
#include <map>

#include <iostream>

namespace quadrature
{
	// Cached storage.
	namespace
	{
		std::map<std::pair<int, int>, double> integrationPoints;
		std::map<std::pair<int, int>, double> weights;
	}

	/******************************************************************************
	 * legendrePolynomial
	 * 
	 * @details    Returns the ith derivative of the nth Legendre polynomial.
	 *
	 * @param[in] n 			Gives the degree of the polynomial.
	 * @param[in] i 			The ith derivative.
	 * @param[in] x 			The point at which to evaluate the polynomial.
	 * @return 					The Legendre polynomial.
	 ******************************************************************************/
	f_double legendrePolynomial(const int &a_n, const int &a_i)
	{
		if (a_i==0)
		{
			switch(a_n)
			{
				case 0: return[](double x)->double
					{
						return 1;
					};
					break;
				case 1: return[](double x)->double
					{
						return x;
					};
					break;
				default: return[a_n](double x)->double
					{
						using namespace common;

						return 	addFunction(
									constantMultiplyFunction(
										double(2*a_n-1)/a_n * x,
										legendrePolynomial(a_n-1, 0)
									),
									constantMultiplyFunction(
										-1,
										constantMultiplyFunction(
											double(a_n-1)/a_n,
											legendrePolynomial(a_n-2, 0)
										)
									)
								)(x);
					};
			}
		}
		else
		{
			switch(a_n)
			{
				case 0: return[](double x)->double
					{
						return 0;
					};
					break;
				case 1: return[a_i](double x)->double
					{
						return (a_i==1)?1:0;
					};
					break;
				default: return[a_n, a_i](double x)->double
					{
						using namespace common;

						return 	addFunction(
									addFunction(
										constantMultiplyFunction(
											double(2*a_n-1)/(a_n-1) * x,
											legendrePolynomial(a_n-1, a_i)
										),
										constantMultiplyFunction(
											-1,
											constantMultiplyFunction(
												double(a_n)/(a_n-1),
												legendrePolynomial(a_n-2, a_i)
											)
										)
									),
									constantMultiplyFunction(
										(a_i-1) * double(2*a_n-1)/(a_n-1),
										legendrePolynomial(a_n-1, a_i-1)
									)
								)(x);
					};
			}
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

		while (fabs(legendrePolynomial(n, 0)(root)) >= 1e-5)
			root = root - legendrePolynomial(n, 0)(root)/legendrePolynomial(n, 1)(root);

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
	double gaussLegendreQuadrature(const f_double f, const int n)
	{
		double weight, x;
		double h = double(2)/n;
		double quadrature = 0;

		for (int i=0; i<n; ++i)
		{
			x = legendrePolynomialRoot(n, i);
			weight = double(2) / ((1 - pow(x, 2))*pow(legendrePolynomial(n, 1)(x), 2));
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

	double get_gaussLegendrePoint(const int &a_n, const int &a_i)
	{
		if (integrationPoints.find(std::make_pair(a_n, a_i)) == integrationPoints.end())
		{
			integrationPoints[std::make_pair(a_n, a_i)] = legendrePolynomialRoot(a_n, a_i);
		}

		return integrationPoints[std::make_pair(a_n, a_i)];
	}

	double get_gaussLegendreWeight(const int &a_n, const int &a_i)
	{
		if (weights.find(std::make_pair(a_n, a_i)) == weights.end())
		{
			double xi = get_gaussLegendrePoint(a_n, a_i);
			weights[std::make_pair(a_n, a_i)] = double(2) / ((1 - pow(xi, 2))*pow(legendrePolynomial(a_n, 1)(xi), 2));
		}

		return weights[std::make_pair(a_n, a_i)];
	}
}