/******************************************************************************
 * @details This is a file containing declarations of the [quadrature]
 * 				namespace.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/12/03
 ******************************************************************************/
#ifndef NAMESPACE_QUADRATURE
#define NAMESPACE_QUADRATURE

#include "common.hpp"
#include <cassert>
#include <cmath>
#include <functional>

namespace quadrature
{
	double   gaussLegendreQuadrature(const f_double f, const int n);
	f_double legendrePolynomial(const int &a_n, const int &a_i);
	double   legendrePolynomialRoot(const int n, const int i);
	void     legendrePolynomialRoots(const int n, double roots[]);
	double   trapeziumRule(const int n, const double fValues[], const double h);

	double get_gaussLegendrePoint(const int &a_n, const int &a_i);
	double get_gaussLegendreWeight(const int &a_n, const int &a_i);
}

#endif