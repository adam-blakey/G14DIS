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
	f_double addFunction(const f_double &f, const f_double &g);
	double** allocateMatrix(const int &noRows, const int &noCols);
	f_double constantMultiplyFunction(const double &a, const f_double &f);
	void     copyArray(const int &n, const double* const &array1, double* const &array2);
	void     deallocateMatrix(const int &noRows, double** &matrix);
	f_double multiplyFunction(const f_double &f, const f_double &g);	
	double   referenceL2Norm(const int &n, const f_double &u);
	f_double transformFunction(const f_double &f, const double &xj, const double &xjp1);
	void     tridiagonalVectorMultiplication(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &c, const std::vector<double> &x, std::vector<double> &solution);
	void     setToZero(const int &n, double* const array);
}

#endif