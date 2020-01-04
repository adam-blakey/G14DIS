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
Matrix_sparse<T>::Matrix_sparse(const int &N)
: Matrix(N, N)
{
	//
}

/******************************************************************************
 * __Matrix__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_sparse<T>::Matrix_sparse(const int &a_noColumns, const int &a_noRows)
{
	this->noColumns = a_noColumns;
	this->noRows = a_noRows;
	items.reserve(noRows * noColumns);
}

/******************************************************************************
 * __Matrix__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_sparse<T>::Matrix_sparse(const int &a_noColumns, const int &a_noRows, const T &a_initial)
: Matrix_sparse(a_noColumns, a_noRows)
{
	for (int i=0; i<noColumns; ++i)
		for (int j=0; j<noRows; ++j)
			items[get_index(i, j)] = a_initial;
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
	for (int i=0; i<noColumns; ++i)
		for (int j=0; j<noRows; ++j)
			items[get_index(i, j)] = a_matrix(i, j);
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

#endif