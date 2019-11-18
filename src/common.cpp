/**
 * @details This is a file containing functions regarding common functions.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/11/11
 */
#include <cmath>
#include <functional>

namespace common
{
	using namespace std;

	/**
	 * __setToZero__
	 * 
	 * @details 	This sets all n elements of the double array to zero.
	 * 
	 * @param[in] n 		The number of elements in the array.
	 * @param[in] array 	A pointer to the first element of the array.
	 */
	void setToZero(const int n, double array[])
	{
		for (int i=0; i<n; ++i)
		{
			array[i] = 0;
		}
	}

	/**
	 * __transformFunction__
	 * 
	 * @details 	Transforms the given function to the reference element.
	 * 
	 * @param[in] f 		The function.
	 * @param[in] n 		Left node of the element.
	 * @param[in] array 	Right node of the element.
	 * @return     			A function that is transformed to the reference element using the left and right side of a 1D element.
	 */
	function<double(double)> transformFunction(function<double(double)> f, double xj, double xjp1)
	{
		return [=](double x)->double{ return f((x+1)*(xjp1-xj)/2 + xj); };
	}

	/**
	 * __addFunction__
	 * 
	 * @details 	Returns the sum of two functions.
	 * 
	 * @param[in] f 		First function
	 * @param[in] g 		Second function.
	 * @return     			A function that is `h(x) = f(x) + g(x)`.
	 */
	function<double(double)> addFunction(function<double(double)> f, function<double(double)> g)
	{
		return [=](double x)->double{ return f(x) + g(x); };
	}

	/**
	 * __multiplyFunction__
	 * 
	 * @details 	Returns the product of two functions.
	 * 
	 * @param[in] f 		First function
	 * @param[in] g 		Second function.
	 * @return     			A function that is `h(x) = f(x) * g(x)`.
	 */
	function<double(double)> multiplyFunction(function<double(double)> f, function<double(double)> g)
	{
		return [=](double x)->double{ return f(x) * g(x); };
	}

	/**
	 * __constantMultiplyFunction__
	 * 
	 * @details 	Returns a the value of a function multiplied by a constant.
	 * 
	 * @param[in] a 		The constant.
	 * @param[in] f 		The function.
	 * @return     			A function that is `h(x) = a * f(x)`.
	 */
	function<double(double)> constantMultiplyFunction(double a, function<double(double)> f)
	{
		return [=](double x)->double{ return a * f(x); };
	}

	/**
	 * __referenceL2Norm__
	 * 
	 * @details 	Returns the L2 norm on the reference element.
	 * 
	 * @param[in] f 		First function
	 * @param[in] g 		Second function.
	 * @return     			A function that is `h(x) = f(x) + g(x)`.
	 */
	double referenceL2Norm(const int n, std::function<double(double)> u)
	{
		double norm = 0;
		double h = double(2)/(n-1);
		double node;

		for (int i=0; i<n; ++i)
		{
			node = -1 + i*h;
			norm += pow(u(node), 2)*h;
		}

		return sqrt(norm);
	}

	/**
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
	 */
	void tridiagonalVectorMultiplication(const int n, const double a[], const double b[], const double c[], const double x[], double solution[])
	{
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
}