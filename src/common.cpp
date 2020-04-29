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

#include <iostream>

namespace common
{
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
	f_double transformFunction(const f_double &a_f, const double &a_xj, const double &a_xjp1)
	{
		return [=](double x)->double{ return a_f((x+1)*(a_xjp1-a_xj)/2 + a_xj); };
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
	f_double addFunction(const f_double &a_f, const f_double &a_g)
	{
		return [=](double x)->double{ return a_f(x) + a_g(x); };
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
	f_double multiplyFunction(const f_double &a_f, const f_double &a_g)
	{
		return [=](double x)->double{ return a_f(x) * a_g(x); };
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
	f_double constantMultiplyFunction(const double &a_a, const f_double &a_f)
	{
		return [=](double x)->double{ return a_a * a_f(x); };
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
	void tridiagonalVectorMultiplication(const std::vector<double> &a_a, const std::vector<double> &a_b, const std::vector<double> &a_c, const std::vector<double> &a_x, std::vector<double> &a_solution)
	{
		int n = a_b.size();
		
		for (int i=0; i<n; ++i)
		{
			if (i==0)
				a_solution[i] = a_b[0]*a_x[0] + a_c[0]*a_x[1];
			else if (i==n-1)
				a_solution[i] = a_a[n-2]*a_x[n-2] + a_b[n-1]*a_x[n-1];
			else
				a_solution[i] = a_a[i-1]*a_x[i-1] + a_b[i]*a_x[i] + a_c[i]*a_x[i+1];
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
		f_double f2 = common::multiplyFunction(a_f, a_f);

		// bad Adam
		return sqrt(quadrature::gaussLegendreQuadrature(f2, 8));
	}

	/******************************************************************************
	 * __L2NormDifference__
	 * 
	 * @details     Calculates the L2 norm of the difference between two functions.
	 *
	 * @param[in]  a_f 	An input function 1.
	 * @param[in]  a_g  An input function 2.
	 *
	 * @return     The value of the L2 norm on the domain [-1, 1].
	 ******************************************************************************/
	double L2NormDifference(const f_double &a_f, const f_double &a_g)
	{
		f_double gm    = common::constantMultiplyFunction(-1, a_g);
		f_double diff  = common::addFunction(a_f, gm);

		return L2Norm(diff);
	}

	/******************************************************************************
	 * __vectorL2Norm__
	 * 
	 * @details     Calculates the vector L2 norm of the function.
	 *
	 * @param[in]  a_v1 Vector 1.
	 * @param[in]  a_v2 Vector 2.
	 *
	 * @return     The value of the vector L2 norm.
	 ******************************************************************************/
	double l2Norm(const std::vector<double> &a_v1, const std::vector<double> &a_v2)
	{
		double norm = 0;

		for (int i=0; i<a_v1.size(); ++i)
		{
			//std::cout << a_v1[i] - a_v2[i] << std::endl;
			norm += pow(a_v1[i] - a_v2[i], 2);
		}

		return sqrt(norm);
	}
}

// OPERATOR-
std::vector<double> operator-(const std::vector<double> &a_v1, const std::vector<double> &a_v2)
{
	std::vector<double> result = a_v1;

	for (int i=0; i<result.size(); ++i)
		result[i] -= a_v2[i];

	return result;
}

// OPERATOR*
std::vector<double> operator*(const double &a_c, const std::vector<double> &a_v)
{
	std::vector<double> result = a_v;

	for (int i=0; i<result.size(); ++i)
		result[i] = a_c*a_v[i];

	return result;
}