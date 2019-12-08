/******************************************************************************
 * @details This is a file containing functions regarding common functions.
 * 
 * @author     Adam Matthew Blakey
 * @date       2019/11/11
 ******************************************************************************/
#include <cmath>
#include <functional>

namespace common
{
	/******************************************************************************
	 * __setToZero__
	 * 
	 * @details 	This sets all n elements of the double array to zero.
	 * 
	 * @param[in] n 		The number of elements in the array.
	 * @param[in] array 	A pointer to the first element of the array.
	 ******************************************************************************/
	void setToZero(const int &n, double* const array)
	{
		for (int i=0; i<n; ++i)
			array[i] = 0;
	}

	/******************************************************************************
	 * __copyArray__
	 * 
	 * @details 	Does a deep copy on an array.
	 * 
	 * @param[in] n 		The number of elements in the arrays.
	 * @param[in] array1 	The source array.
	 * @param[in] array2 	The destination array.
	 ******************************************************************************/
	void setToZero(const int &n, const double* const &array1, double* const &array2)
	{
		for (int i=0; i<n; ++i)
			array2[i] = array1[i];
	}

	/******************************************************************************
	 * __allocateMatrix__
	 * 
	 * @details 	Allocates a matrix given the number of rows and columns.
	 * 
	 * @param[in] noRows 	Number of rows.
	 * @param[in] noCols 	Number of columns.
	 * @return  			A double double-pointer.
	 ******************************************************************************/
	double** allocateMatrix(const int &noRows, const int &noCols)
	{
		double** matrix = new double*[noRows];
		for (int i=0; i<noRows; ++i)
			matrix[i] = new double[noCols];

		return matrix;
	}

	/******************************************************************************
	 * __deallocateMatrix__
	 * 
	 * @details 	Deallocates a matrix, given the number of rows in that matrix.
	 * 
	 * @param[in] noRows 	Number of rows.
	 * @param[out] matrix 	Matrix to deallocate.
	 ******************************************************************************/
	void deallocateMatrix(const int &noRows, double** &matrix)
	{
		for (int i=0; i<noRows; ++i)
			delete[] matrix[i];

		delete[] matrix;
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
	std::function<double(double)> transformFunction(const std::function<double(double)> &f, const double &xj, const double &xjp1)
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
	std::function<double(double)> addFunction(const std::function<double(double)> &f, const std::function<double(double)> &g)
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
	std::function<double(double)> multiplyFunction(const std::function<double(double)> &f, const std::function<double(double)> &g)
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
	std::function<double(double)> constantMultiplyFunction(const double &a, const std::function<double(double)> &f)
	{
		return [=](double x)->double{ return a * f(x); };
	}

	/******************************************************************************
	 * __referenceL2Norm__
	 * 
	 * @details 	Returns the L2 norm on the reference element.
	 * 
	 * @param[in] f 		First function
	 * @param[in] g 		Second function.
	 * @return     			A function that is `h(x) = f(x) + g(x)`.
	 ******************************************************************************/
	double referenceL2Norm(const int &n, const std::function<double(double)> &u)
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
	void tridiagonalVectorMultiplication(const int &n, const double* const &a, const double* const &b, const double* const &c, const double* const &x, double* const &solution)
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