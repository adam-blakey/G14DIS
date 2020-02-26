/******************************************************************************
 * @details This is a file containing declarations of [Matrix].
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/01/02
 ******************************************************************************/
#ifndef CLASS_SRC_MATRIX_FULL
#define CLASS_SRC_MATRIX_FULL

#include "common.hpp"
#include "matrix.hpp"
#include "matrix_full.hpp"
#include <cmath>
#include <iostream>

/******************************************************************************
 * __Matrix_full__
 * 
 * @details 	The default [Matrix] constructor.
 ******************************************************************************/
template<class T>
Matrix_full<T>::Matrix_full(const int &N)
: Matrix_full(N, N)
{
	//
}

/******************************************************************************
 * __Matrix_full__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_full<T>::Matrix_full(const int &a_noColumns, const int &a_noRows)
{
	this->noColumns = a_noColumns;
	this->noRows = a_noRows;
	items.resize(this->noRows * this->noColumns);
}

/******************************************************************************
 * __Matrix_full__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_full<T>::Matrix_full(const int &a_noColumns, const int &a_noRows, const T &a_initial)
: Matrix_full(a_noColumns, a_noRows)
{
	for (int i=0; i<this->noColumns; ++i)
		for (int j=0; j<this->noRows; ++j)
			items[get_index(i, j)] = a_initial;
}

/******************************************************************************
 * __Matrix_full__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_full<T>::Matrix_full(const Matrix<T> &a_matrix)
: Matrix_full(a_matrix.get_noColumns(), a_matrix.get_noRows())
{
	for (int i=0; i<this->noColumns; ++i)
		for (int j=0; j<this->noRows; ++j)
			items[get_index(i, j)] = a_matrix(i, j);
}

/******************************************************************************
 * __resize__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
void Matrix_full<T>::resize(const int &a_noNonZeros)
{
	items.resize(a_noNonZeros);
}

/******************************************************************************
 * __get_index__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
int Matrix_full<T>::get_index(const int &a_x, const int &a_y) const
{
	if (a_x >= this->noColumns || a_y >= this->noRows)
	{
		std::cerr << "Error: Requested indices exceed matrix dimensions." << std::endl;
		return 0;
	}

	return a_x + a_y*this->noColumns;
}

/******************************************************************************
 * __item__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
T& Matrix_full<T>::item(const int &a_x, const int &a_y)
{
	// Dimensions must be the same.
	if (a_x >= this->noColumns || a_y >= this->noRows)
	{
		std::cerr << "Invalid index at (" << a_x << ", " << a_y << ").";
		return items[0];
	}

	return items[get_index(a_x, a_y)];
}

/******************************************************************************
 * __item__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
const T Matrix_full<T>::item(const int &a_x, const int &a_y) const
{
	// Dimensions must be the same.
	if (a_x >= this->noColumns || a_y >= this->noRows)
	{
		std::cerr << "Invalid index.";
		return items[0];
	}

	return items[get_index(a_x, a_y)];
}

/******************************************************************************
 * __operator+__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_full<T> Matrix_full<T>::operator+(const Matrix<T> &a_RHS)
{
	// Creates new matrix and calculates elements appropriately.
	Matrix_full<T> tempMatrix((*this));
	
	tempMatrix += a_RHS;

	return tempMatrix;
}

/******************************************************************************
 * __operator-__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_full<T> Matrix_full<T>::operator-(const Matrix<T> &a_RHS)
{
	// Creates new matrix and calculates elements appropriately.
	Matrix_full<T> tempMatrix((*this));
	
	tempMatrix -= a_RHS;

	return tempMatrix;
}

/******************************************************************************
 * __operator*__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_full<T> Matrix_full<T>::operator*(const Matrix<T> &a_RHS)
{
		// Matching dimensions.
	if (this->noRows != a_RHS.get_noColumns())
	{
		std::cerr << "Matrix dimensions do not match (cannot multiply "
			<< this->noColumns
			<< "x"
			<< this->noRows
			<< " by "
			<< a_RHS.get_noColumns()
			<< "x"
			<< a_RHS.get_noRows()
			<< ").";

		return *this;
	}

	int newRows = a_RHS.get_noRows();
	int newColumns = this->noColumns;

	Matrix_full<T> tempMatrix(newColumns, newRows, 0);

	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<newColumns; ++i)
		for (int j=0; j<newRows; ++j)
			for (int k=0; k<this->get_noRows(); ++k)
				tempMatrix(i, j) += item(i, k) * a_RHS(k, j);

	return tempMatrix;
}

/******************************************************************************
 * __operator*__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_full<T> Matrix_full<T>::operator*(const T &a_RHS)
{
	Matrix_full<T> tempMatrix(*this);

	tempMatrix *= a_RHS;

	return tempMatrix;
}

/******************************************************************************
 * __operator/__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_full<T> Matrix_full<T>::operator/(const T &a_RHS)
{
	Matrix_full<T> tempMatrix(*this);

	tempMatrix /= a_RHS;

	return tempMatrix;
}

/******************************************************************************
 * __get_noRows__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
int Matrix_full<T>::get_noRows() const
{
	return this->noRows;
}

/******************************************************************************
 * __get_noColumns__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
int Matrix_full<T>::get_noColumns() const
{
	return this->noColumns;
}

/******************************************************************************
 * __set__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
void Matrix_full<T>::set(const int &a_x, const int &a_y, const T &a_value)
{
	item(a_x, a_y) = a_value;
}

/******************************************************************************
 * STANDALONE VECTOR TENSORISATION -- IS THIS THE BEST PLACE FOR IT?!?!?!
******************************************************************************/
template<class T>
Matrix_full<T> Tensorise(const std::vector<T> &a_vec1, const std::vector<T> &a_vec2)
{
	// Temporary matrices.
	Matrix_full<T> mat1(a_vec1.size(), 1);
	for (int i=0; i<a_vec1.size(); ++i)
		mat1(1, i) = a_vec1[i];

	Matrix_full<T> mat2(1, a_vec2.size());
	for (int i=0; i<a_vec2.size(); ++i)
		mat2(i, 1) = a_vec2[i];

	// Matrix result.
	Matrix_full<T> result = mat1 * mat2;

	return result;
}

/******************************************************************************
 * STANDALONE VECTOR REDUCTION -- IS THIS THE BEST PLACE FOR IT?!?!?!
******************************************************************************/
template<class T>
T Scalarise(const std::vector<T> &a_vec1, const std::vector<T> &a_vec2)
{
	// Temporary matrices.
	Matrix_full<T> mat1(1, a_vec1.size());
	for (int i=0; i<a_vec1.size(); ++i)
		mat1(1, i) = a_vec1[i];

	Matrix_full<T> mat2(a_vec2.size(), 1);
	for (int i=0; i<a_vec2.size(); ++i)
		mat2(i, 1) = a_vec2[i];

	// Single value result.
	Matrix_full<T> result = mat1 * mat2;

	return result(0, 0);
}

/******************************************************************************
 * STANDALONE VECTOR NORM -- IS THIS THE BEST PLACE FOR IT?!?!?!
******************************************************************************/
template<class T>
T l2Norm(const std::vector<T> &a_vec)
{
	T value = 0;

	for (int i=0; i<a_vec.size(); ++i)
		value += pow(a_vec[i], 2);

	return sqrt(value);
}

/******************************************************************************
 * __innerProduct__
 * 
 * @details 	Calculates the inner product of two vectors.
 *
 * @param[in]  a_vec1  First vector.
 * @param[in]  a_vec2  Second vector.
 *
 * @tparam     General type.
 *
 * @return     Returns the sum of the product of individual vector components.
 ******************************************************************************/
template<class T>
T innerProduct(const std::vector<T> &a_vec1, const std::vector<T> &a_vec2)
{
	T result = 0;

	for (int i=0; i<a_vec1.size(); ++i)
		result += a_vec1[i] * a_vec2[i];

	return result;
}

/******************************************************************************
 * __linearCombination__
 * 
 * @details 	Calculates a linear combination of two vectors.
 *
 * @param[in]  a_vec1  First vector.
 * @param[in]  a_vec2  Second vector.
 *
 * @tparam     General type.
 *
 * @return     Returns the linear combination.
 ******************************************************************************/
template<class T>
std::vector<T> linearCombination(const T &a_constant1, const std::vector<T> &a_vec1, const T &a_constant2, const std::vector<T> &a_vec2)
{
	std::vector<T> result(a_vec1.size());

	for (int i=0; i<result.size(); ++i)
		result[i] = a_constant1*a_vec1[i] + a_constant2*a_vec2[i];

	return result;
}

/******************************************************************************
 * STANDALONE CONJUGATE GRADIENT -- IS THIS THE BEST PLACE FOR IT?!?!?!
******************************************************************************/
template<class T>
std::vector<T> ConjugateGradient(const Matrix<T> &A, const std::vector<T> &b)
{
	const double NEARZERO = 1e-10;
	const double TOLERANCE = 1e-10;

   	int n = A.get_noRows();
	std::vector<T> X(n, 0);

	std::vector<T> R = b;
	std::vector<T> P = R;
	int k = 0;

	while (k < n)
	{
	  std::vector<T> Rold = R;
	  std::vector<T> AP = A*P;

	  double alpha = innerProduct(R, R)/ std::max(innerProduct(P, AP), NEARZERO);
	  X = linearCombination(1.0, X, alpha, P);            // Next estimate of solution
	  R = linearCombination(1.0, R, -alpha, AP);          // Residual 

	  if (l2Norm(R) < TOLERANCE) break;             // Convergence test

	  double beta = innerProduct(R, R) / std::max(innerProduct(Rold, Rold), NEARZERO);
	  P = linearCombination(1.0, R, beta, P);             // Next gradient

	  ++k;
	}

	return X;
}


/*template<class T>
std::vector<T> ConjugateGradient(const Matrix<T> &a_matrix, const std::vector<T> &a_vector)
{
	// Tolerance.
	const int TOL = 1e-10;

	// Sets up some temporary variables.
	std::vector<T> y_m  (a_matrix.get_noColumns());
	std::vector<T> y_mm1(a_matrix.get_noColumns(), 0); // I.C.
	std::vector<T> r_m  (a_matrix.get_noColumns());
	std::vector<T> r_mm1(a_matrix.get_noColumns());
	std::vector<T> p_m  (a_matrix.get_noColumns());
	std::vector<T> p_mp1(a_matrix.get_noColumns());

	// Initialises the temporary variables.
	r_mm1 = a_vector - (a_matrix*y_mm1);
	p_m = r_mm1;

	// Loops while the norm is big.
	while(l2Norm(r_mm1) > TOL)
	{
		// Applies conjugate gradient.
		double alpha = Scalarise(r_mm1, r_mm1)/Scalarise(p_m, a_matrix*p_m);
		y_m = y_mm1 + alpha*p_m;
		r_m = r_mm1 - alpha*a_matrix*p_m;
		double beta = Scalarise(r_m, r_m)/Scalarise(r_mm1, r_mm1);
		p_mp1 = r_m + beta*p_m;

		// Sets variables for next iteration.
		r_mm1 = r_m;
		y_mm1 = y_m;
		p_m = p_mp1;
	}

	return y_m;
}*/

#endif