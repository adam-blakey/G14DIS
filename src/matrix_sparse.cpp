/******************************************************************************
 * @details This is a file containing declarations of [Matrix].
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/01/02
 ******************************************************************************/
#ifndef CLASS_SRC_MATRIX_SPARSE
#define CLASS_SRC_MATRIX_SPARSE

#include "matrix.hpp"
#include "matrix_sparse.hpp"
#include <cmath>
#include <iostream>

/******************************************************************************
 * __Matrix__
 * 
 * @details 	The default [Matrix] constructor.
 ******************************************************************************/
template<class T>
Matrix_sparse<T>::Matrix_sparse(const int &a_N, const int &a_noNonZero)
: Matrix(a_N, a_N, a_noNonZero)
{
	//
}

/******************************************************************************
 * __Matrix__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_sparse<T>::Matrix_sparse(const int &a_noColumns, const int &a_noRows, const int &a_noNonZero)
{
	this->noColumns = a_noColumns;
	this->noRows = a_noRows;

	this->resize(a_noNonZero);
}

/******************************************************************************
 * __Matrix__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_sparse<T>::Matrix_sparse(const Matrix<T> &a_matrix)
: Matrix_sparse(a_matrix.get_noColumns(), a_matrix.get_noRows())
{
	// Temporary vectors with matrix entries and column numbers.
	std::vector<T> tempMatrixEntries(0);
	std::vector<T> tempColumnNos(0);

	// Resizes the row starts vector to correct length.
	this->rowStarts.resize(this->noRows + 1);

	for (int j=0; j<noRows; ++j)
	{
		bool rowStartAssigned = false;

		for (int i=0; i<noColumns; ++i)
		{
			// Insert element if nonzero.
			if (a_matrix(i, j) != 0)
			{
				tempMatrixEntries.push_back(a_matrix(i, j));
				tempColumnNos    .push_back(i);

				if (!rowStartAssigned)
				{
					this->rowStarts[j] = tempMatrixEntries.size();
					rowStartAssigned = true;
				}
			}
		}

		// If the row start still hasn't been assigned.
		if (!rowStartAssigned)
			rowStarts[j] = tempMatrixEntries.size();
	}

	this->rowStarts[noRows] = tempColumnsNos.size();

	this->matrixEntries = tempMatrixEntries;
	this->columnNos = tempColumnNos;
}

/******************************************************************************
 * __resize__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
void Matrix_sparse<T>::resize(const int &a_noNonZeros) const
{
	matrixEntries.resize(a_noNonZeros);
	columnNos    .resize(a_noNonZeros);
	rowStarts    .resize(noRows+1);
}

/******************************************************************************
 * __get_index__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
int Matrix_sparse<T>::get_index(const int &a_x, const int &a_y) const
{
	if (a_x >= noColumns || a_y >= noRows)
	{
		std::cerr << "Error: Requested indices exceed matrix dimensions." << std::endl;
		return 0;
	}

	return a_x + a_y*noColumns;
}

/******************************************************************************
 * __item__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
T& Matrix_sparse<T>::item(const int &a_x, const int &a_y)
{
	// Dimensions must be the same.
	if (a_x >= noColumns || a_y >= noRows)
	{
		std::cerr << "Invalid index.";
		return items[8];
	}

	return items[get_index(a_x, a_y)];
}

/******************************************************************************
 * __item__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
const T& Matrix_sparse<T>::item(const int &a_x, const int &a_y) const
{
	// Dimensions must be the same.
	if (a_x >= noColumns || a_y >= noRows)
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
Matrix_sparse<T> Matrix_sparse<T>::operator+(const Matrix<T> &a_RHS)
{
	// Creates new matrix and calculates elements appropriately.
	Matrix_sparse<T> tempMatrix((*this));
	
	tempMatrix += a_RHS;

	return tempMatrix;
}

/******************************************************************************
 * __operator-__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_sparse<T> Matrix_sparse<T>::operator-(const Matrix<T> &a_RHS)
{
	// Creates new matrix and calculates elements appropriately.
	Matrix_sparse<T> tempMatrix((*this));
	
	tempMatrix -= a_RHS;

	return tempMatrix;
}

/******************************************************************************
 * __operator*__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_sparse<T> Matrix_sparse<T>::operator*(const Matrix<T> &a_RHS)
{
		// Matching dimensions.
	if (noRows != a_RHS.get_noColumns())
	{
		std::cerr << "Matrix dimensions do not match (cannot multiply "
			<< noColumns
			<< "x"
			<< noRows
			<< " by "
			<< a_RHS.get_noColumns()
			<< "x"
			<< a_RHS.get_noRows()
			<< ").";

		return *this;
	}

	int newRows = a_RHS.get_noRows();
	int newColumns = noColumns;

	Matrix_sparse<T> tempMatrix(newColumns, newRows, 0);

	// Creates new matrix and calculates elements appropriately.
	for (int i=0; i<newColumns; ++i)
		for (int j=0; j<newRows; ++j)
			for (int k=0; k<noRows; ++k)
				tempMatrix(i, j) += item(i, k) * a_RHS(k, j);

	return tempMatrix;
}

/******************************************************************************
 * __operator*__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_sparse<T> Matrix_sparse<T>::operator*(const T &a_RHS)
{
	Matrix_sparse<T> tempMatrix(*this);

	tempMatrix *= a_RHS;

	return tempMatrix;
}

/******************************************************************************
 * __operator/__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_sparse<T> Matrix_sparse<T>::operator/(const T &a_RHS)
{
	Matrix_sparse<T> tempMatrix(*this);

	tempMatrix /= a_RHS;

	return tempMatrix;
}

#endif