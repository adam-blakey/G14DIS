/******************************************************************************
 * @details This is a file containing declarations of [Matrix].
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/01/02
 ******************************************************************************/
#ifndef CLASS_SRC_MATRIX_FULL
#define CLASS_SRC_MATRIX_FULL

#include "matrix.hpp"
#include "matrix_full.hpp"
#include <cmath>
#include <iostream>

/******************************************************************************
 * __Matrix__
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
 * __Matrix__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_full<T>::Matrix_full(const int &a_noColumns, const int &a_noRows)
{
	this->noColumns = a_noColumns;
	this->noRows = a_noRows;
	items.reserve(this->noRows * this->noColumns);
}

/******************************************************************************
 * __Matrix__
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
 * __Matrix__
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
const T& Matrix_full<T>::item(const int &a_x, const int &a_y) const
{
	// Dimensions must be the same.
	if (a_x >= this->noColumns || a_y >= this->noRows)
	{
		std::cerr << "Invalid index.";
		return items[0];
	}

	return items[get_index(a_x, a_y)];
}

#endif