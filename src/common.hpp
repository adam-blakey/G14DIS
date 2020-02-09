/******************************************************************************
 * @details Declarations for [Common] namespace.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/11/11
 ******************************************************************************/
#ifndef NAMESPACE_COMMON
#define NAMESPACE_COMMON

#include <cmath>
#include <functional>
#include <vector>

typedef std::function<double(double)> f_double;

namespace common
{
	f_double addFunction(const f_double &a_f, const f_double &a_g);
	f_double constantMultiplyFunction(const double &a_a, const f_double &a_f);
	double   L2Norm(const f_double &a_f);
	double   L2NormDifference(const f_double &a_f, const f_double &a_g);
	f_double multiplyFunction(const f_double &a_f, const f_double &a_g);	
	f_double transformFunction(const f_double &a_f, const double &a_xj, const double &a_xjp1);
	void     tridiagonalVectorMultiplication(const std::vector<double> &a_a, const std::vector<double> &a_b, const std::vector<double> &a_c, const std::vector<double> &a_x, std::vector<double> &a_solution);
}

#endif