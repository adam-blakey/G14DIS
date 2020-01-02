#ifndef NAMESPACE_QUADRATURE
#define NAMESPACE_QUADRATURE

#include "common.hpp"
#include <cassert>
#include <cmath>
#include <functional>

namespace quadrature
{
	f_double legendrePolynomial(const int n);

	f_double legendrePolynomialDerivative(const int n);

	double legendrePolynomialRoot(const int n, const int i);

	void legendrePolynomialRoots(const int n, double roots[]);
	
	double gaussLegendreQuadrature(const f_double f, const int n);
	
	double trapeziumRule(const int n, const double fValues[], const double h);
}

#endif