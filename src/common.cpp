/******************************************************************************
 * @details This is a file containing functions regarding [common] functions.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/11/11
 ******************************************************************************/
#include "common.hpp"
#include "quadrature.hpp"
#include <cmath>
#include <functional>

namespace common
{
	/******************************************************************************
	 * __copyArray__
	 * 
	 * @details 	Does a deep copy on an array.
	 * 
	 * @param[in] n 		The number of elements in the arrays.
	 * @param[in] array1 	The source array.
	 * @param[in] array2 	The destination array.
	 ******************************************************************************/
	void copyArray(const int &n, const double* const &array1, double* const &array2)
	{
		for (int i=0; i<n; ++i)
			array2[i] = array1[i];
	}

	/******************************************************************************
	 * __transformFunction__
	 * 
	 * @details 	Transforms the given function to the reference element.
	 * 
	 * @param[in] f 		The function.
	 * @param[in] xj 		The left side of the element.
	 * @param[in] xjp1 		The right side of the element.
	 * @return     			A function that is transformed to the reference element using the left and right side of a 1D element.
	 ******************************************************************************/
	f_double transformFunction(const f_double &f, const double &xj, const double &xjp1)
	{
		return [=](double x)->double{ return f((x+1)*(xjp1-xj)/2 + xj); };
	}

	/******************************************************************************
	 * __addFunction__
	 * 
	 * @details 	Returns the sum of two functions.
	 * 
	 * @param[in] f 		First function
	 * @param[in] g 		Second function.
	 * @return     			A function that is `h(x) = f(x) + g(x)`.
	 ******************************************************************************/
	f_double addFunction(const f_double &f, const f_double &g)
	{
		return [=](double x)->double{ return f(x) + g(x); };
	}

	/******************************************************************************
	 * __multiplyFunction__
	 * 
	 * @details 	Returns the product of two functions.
	 * 
	 * @param[in] f 		First function
	 * @param[in] g 		Second function.
	 * @return     			A function that is `h(x) = f(x) * g(x)`.
	 ******************************************************************************/
	f_double multiplyFunction(const f_double &f, const f_double &g)
	{
		return [=](double x)->double{ return f(x) * g(x); };
	}

	/******************************************************************************
	 * __constantMultiplyFunction__
	 * 
	 * @details 	Returns a the value of a function multiplied by a constant.
	 * 
	 * @param[in] a 		The constant.
	 * @param[in] f 		The function.
	 * @return     			A function that is `h(x) = a * f(x)`.
	 ******************************************************************************/
	f_double constantMultiplyFunction(const double &a, const f_double &f)
	{
		return [=](double x)->double{ return a * f(x); };
	}

	/******************************************************************************
	 * __tridiagonalVectorMultiplication__
	 * 
	 * @details 	Calculates the product `Ax`.
	 * 
	 * @param[in] n 		The dimension of the square matrix, `A`.
	 * @param[in] a 		Lower diagonal of `A`.
	 * @param[in] b 		Diagonal of `A`.
	 * @param[in] c 		Upper diagonal of `A`.
	 * @param[in] x 		The vector to multiply by.
	 * @param[out] solution The result of the operation.
	 ******************************************************************************/
	void tridiagonalVectorMultiplication(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &c, const std::vector<double> &x, std::vector<double> &solution)
	{
		int n = b.size();
		
		for (int i=0; i<n; ++i)
		{
			if (i==0)
				solution[i] = b[0]*x[0] + c[0]*x[1];
			else if (i==n-1)
				solution[i] = a[n-2]*x[n-2] + b[n-1]*x[n-1];
			else
				solution[i] = a[i-1]*x[i-1] + b[i]*x[i] + c[i]*x[i+1];
		}
	}
	/******************************************************************************
	 * __L2Norm__
	 * 
	 * @details     Calculates the L2 norm of the function.
	 *
	 * @param[in]  a_f 	An input function.
	 *
	 * @return     The value of the L2 norm on the domain [-1, 1].
	 ******************************************************************************/
	double L2Norm(const f_double &a_f)
	{
		quadrature::gaussLegendreQuadrature(a_f, 8);

		return 0;
	}
}