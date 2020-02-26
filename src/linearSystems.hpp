/******************************************************************************
 * @details Declarations for [linearSystems] namespace.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/11/11
 ******************************************************************************/
#ifndef NAMESPACE_LINEARSYSTEMS
#define NAMESPACE_LINEARSYSTEMS

#include "matrix.hpp"
#include <cmath>
#include <functional>
#include <vector>

typedef std::function<double(double)> f_double;

namespace linearSystems
{
	void thomasInvert(const std::vector<double> a, const std::vector<double> b, const std::vector<double> c, const std::vector<double> d, std::vector<double> &solution);
	std::vector<double> conjugateGradient(const Matrix<double> &a_matrix, const std::vector<double> &a_b, const double &a_tolerance);
	double dotProduct(const std::vector<double> &a_v1, const std::vector<double> &a_v2);
}

std::vector<double> operator*(const double &a_constant, const std::vector<double> &a_vector);
std::vector<double>& operator+=(std::vector<double> &a_v1, const std::vector<double> &a_v2);
std::vector<double> operator+(const std::vector<double> &a_v1, const std::vector<double> &a_v2);

#endif