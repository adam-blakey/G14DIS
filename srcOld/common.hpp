/******************************************************************************
 * @details Declarations for [common] namespace.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/11/11
 ******************************************************************************/
#ifndef NAMESPACE_COMMON
#define NAMESPACE_COMMON

#include <cmath>
#include <functional>

typedef std::function<double(double)> f_double;

namespace common
{
	void setToZero(const int &n, double* const array);

	void copyArray(const int &n, const double* const &array1, double* const &array2);
	
	double** allocateMatrix(const int &noRows, const int &noCols);

	void deallocateMatrix(const int &noRows, double** &matrix);
	
	f_double transformFunction(const f_double &f, const double &xj, const double &xjp1);

	f_double addFunction(const f_double &f, const f_double &g);

	f_double multiplyFunction(const f_double &f, const f_double &g);

	f_double constantMultiplyFunction(const double &a, const f_double &f);

	double referenceL2Norm(const int &n, const f_double &u);

	void tridiagonalVectorMultiplication(const int &n, const double* const &a, const double* const &b, const double* const &c, const double* const &x, double* const &solution);
}

#endif