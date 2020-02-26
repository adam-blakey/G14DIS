/******************************************************************************
 * @details This is a file containing declarations of [Matrix].
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/01/02
 ******************************************************************************/
#ifndef CLASS_SRC_MATRIX
#define CLASS_SRC_MATRIX

#include "matrix.hpp"
#include <cmath>
#include <iostream>

/******************************************************************************
 * __get_diagonal__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
std::vector<T> Matrix<T>::get_diagonal() const
{
	int diagonalLength = this->get_noColumns()<this->get_noRows()?this->get_noColumns():this->get_noRows();
	std::vector<T> diagonal(diagonalLength, 0);

	for (int i=0; i<diagonalLength; ++i)
		diagonal[i] = item(i, i);

	return diagonal;
}

/******************************************************************************
 * __operator()__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
const T Matrix<T>::operator()(const int &a_x, const int &a_y) const
{
	return item(a_x, a_y);
}

/******************************************************************************
 * __operator=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &a_RHS)
{
	if (&a_RHS == this)
		return *this;

	this->noRows = a_RHS.get_noRows();
	this->noColumns = a_RHS.get_noColumns();

	this->resize(this->get_noRows() * this->get_noColumns());

	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			item(i, j) = a_RHS(i, j);

	return *this;
}

/******************************************************************************
 * __operator+=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T> &a_RHS)
{
	// Dimensions must be the same.
	if (this->get_noRows() != a_RHS.get_noRows() || this->get_noColumns() != a_RHS.get_noColumns())
	{
		std::cerr << "Matrix dimensions do not match.";
		return *this;
	}

	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			item(i, j) += a_RHS(i, j); // I think it's a problem with this LHS -- I don't think you can write to it...

	return *this;
}

/******************************************************************************
 * __operator-=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T> &a_RHS)
{
	// Dimensions must be the same.
	if (this->get_noRows() != a_RHS.get_noRows() || this->get_noColumns() != a_RHS.get_noColumns())
	{
		std::cerr << "Matrix dimensions do not match.";
		return *this;
	}

	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			item(i, j) -= a_RHS(i, j); // I think it's a problem with this LHS -- I don't think you can write to it...

	return *this;
}

/******************************************************************************
 * __operator*=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T> &a_RHS)
{
	// Matching dimensions.
	if (this->get_noRows() != a_RHS.get_noColumns() || this->get_noRows() != a_RHS.get_noRows())
	{
		std::cerr << "Error: Matrix dimensions do not match." << std::endl;
		return *this;
	}

	Matrix<T> tempMatrix(this->get_noColumns(), this->get_noRows(), 0);

	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			for (int k=0; k<this->get_noRows(); ++k)
				tempMatrix(i, j) += item(i, k) * a_RHS(k, j);

	*this = tempMatrix;

	return *this;
}

/******************************************************************************
 * __operator*=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator*=(const T &a_RHS)
{
	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			item(i, j) *= a_RHS;

	return *this;
}

/******************************************************************************
 * __operator*__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
std::vector<T> Matrix<T>::operator*(const std::vector<T> &a_RHS) const
{
	if (get_noColumns() != a_RHS.size())
	{
		std::cerr << "Error: Matrix-vector dimensions do not match." << std::endl;
		return a_RHS;
	}

	std::vector<T> tempVector(get_noColumns(), 0);
	for (int i=0; i<get_noColumns(); ++i)
		for (int j=0; j<get_noRows(); ++j)
			tempVector[i] += item(i, j) * a_RHS[j];

	return tempVector;
}

template<class T>
std::vector<T> operator*(const Matrix<T> &a_matrix, const std::vector<T> &a_vector)
{
	return a_matrix*a_vector;
}


/******************************************************************************
 * __operator/=__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix<T>& Matrix<T>::operator/=(const T &a_RHS)
{
	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<this->get_noColumns(); ++i)
		for (int j=0; j<this->get_noRows(); ++j)
			item(i, j) /= a_RHS;

	return *this;
}

#endif