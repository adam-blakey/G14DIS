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
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_sparse<T>::Matrix_sparse(const int &a_noNonZero, const int &a_noRows, const int &a_noColumns)
{
	this->resize(a_noNonZero, a_noRows);
	this->noColumns = a_noColumns;
}

/******************************************************************************
 * __Matrix__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
Matrix_sparse<T>::Matrix_sparse(const Matrix<T> &a_matrix)
: Matrix_sparse(a_matrix.get_noNonZero(), a_matrix.get_noRows(), a_matrix.get_noColumns())
{
	// Temporary vectors with matrix entries and column numbers.
	std::vector<T> tempMatrixEntries(0);
	std::vector<T> tempColumnNos(0);

	for (int j=0; j<this->get_noRows(); ++j)
	{
		bool rowStartAssigned = false;

		for (int i=0; i<this->get_noColumns(); ++i)
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

	this->rowStarts[this->get_noRows()] = tempColumnNos.size();

	this->matrixEntries = tempMatrixEntries;
	this->columnNos = tempColumnNos;
}

/******************************************************************************
 * __resize__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
void Matrix_sparse<T>::resize(const int &a_noNonZeros)
{
	this->resize(a_noNonZeros, this->get_noRows());
}

/******************************************************************************
 * __resize__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
void Matrix_sparse<T>::resize(const int &a_noNonZeros, const int &a_noRows)
{
	matrixEntries.resize(a_noNonZeros);
	columnNos    .resize(a_noNonZeros);
	rowStarts    .resize(a_noRows+1);
}

/******************************************************************************
 * __resize__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
void Matrix_sparse<T>::resize(const int &a_noNonZeros, const int &a_noRows, const int &a_noColumns)
{
	this->resize(a_noNonZeros, a_noRows);
	this->noColumns = a_noColumns;
}

/******************************************************************************
 * __get_index__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
int Matrix_sparse<T>::get_index(const int &a_x, const int &a_y) const
{
	// Check indices.
	if (a_x >= this->get_noColumns() || a_y >= this->get_noRows())
	{
		std::cerr << "Error: Requested indices exceed matrix dimensions." << std::endl;
		return 0;
	}

	int thisRowStart = this->rowStarts[a_y];
	int nextRowStart = this->rowStarts[a_y+1];

	int index = thisRowStart;
	while(index<nextRowStart && a_x!=this->columnNos[index])
	{
		++index;
	}

	if (a_x==this->columnNos[index] && thisRowStart!=nextRowStart)
		return index;
	else
		return -1; // Code for a zero.
}

/******************************************************************************
 * __item__
 * 
 * @details 	
 ******************************************************************************/
template<>
inline const int Matrix_sparse<int>::item(const int &a_x, const int &a_y) const
{
	int index = this->get_index(a_x, a_y);

	if (index == -1)
		return 0;
	else
		return matrixEntries[index];
}

/******************************************************************************
 * __item__
 * 
 * @details 	
 ******************************************************************************/
template<>
inline const double Matrix_sparse<double>::item(const int &a_x, const int &a_y) const
{
	int index = this->get_index(a_x, a_y);

	if (index == -1)
		return 0;
	else
		return matrixEntries[index];
}

/******************************************************************************
 * __item__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
const T Matrix_sparse<T>::item(const int &a_x, const int &a_y) const
{
	int index = this->get_index(a_x, a_y);

	if (index == -1)
		return matrixEntries[0]; // This needs to change...
	else
		return matrixEntries[index];
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
	if (this->get_noRows() != a_RHS.get_noColumns())
	{
		std::cerr << "Matrix dimensions do not match (cannot multiply "
			<< this->get_noColumns()
			<< "x"
			<< this->get_noRows()
			<< " by "
			<< a_RHS.get_noColumns()
			<< "x"
			<< a_RHS.get_noRows()
			<< ").";

		return *this;
	}

	int newRows = a_RHS.get_noRows();
	int newColumns = this->get_noColumns();

	Matrix_sparse<T> tempMatrix(newColumns, newRows, 0);

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

/******************************************************************************
 * __get_noNonZero__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
int Matrix_sparse<T>::get_noNonZero() const
{
	return this->matrixEntries.size();
}

/******************************************************************************
 * __get_noRows__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
int Matrix_sparse<T>::get_noRows() const
{
	return this->rowStarts.size() - 1;
}

/******************************************************************************
 * __get_noColumns__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
int Matrix_sparse<T>::get_noColumns() const
{
	return this->noColumns;
}

/******************************************************************************
 * __set__
 * 
 * @details 	
 ******************************************************************************/
template<class T>
void Matrix_sparse<T>::set(const int &a_x, const int &a_y, const T &a_value) // There's a bug somewhere in here...
{
	// If the element already exists, then we just overwrite it.
	if (this->get_index(a_x, a_y) != -1)
	{
		this->matrixEntries[this->get_index(a_x, a_y)] = a_value;
	}
	//this feels bad...
	else if ((this->matrixEntries.size() == 1) and (this->matrixEntries[0] == 0))
	{
		this->matrixEntries[0] = a_value;
		this->columnNos[0] = a_x;

		for (int i=a_y+1; i<this->rowStarts.size(); ++i)
			++this->rowStarts[i];
	}
	else
	{
		int thisRowStart = this->rowStarts[a_y];
		int nextRowStart = this->rowStarts[a_y+1];

		int index;
		bool foundIndex = false;

		for (index=thisRowStart; index<nextRowStart && !foundIndex; ++index)
			if (columnNos[index] < a_x)
				foundIndex = true;

		         std::vector<int>::iterator colIt = this->columnNos    .begin() + index;
		typename std::vector<T>  ::iterator matIt = this->matrixEntries.begin() + index;



		// Inserts value at correct location.
		this->columnNos.insert(colIt, a_x);
		this->matrixEntries.insert(matIt, a_value);

		// Update all following row starts.
		for (int i=a_y+1; i<this->rowStarts.size(); ++i)
			++this->rowStarts[i];

		// DOESN'T QUITE WORK
		// PRINT OUT ALL THREE VECTORS TO DEBUG
	



		// this->matrixEntries
		// this->columnNos
		// this->rowStarts
		 
		/*int rowEntryIndex     = this->rowStarts[a_y];
		int nextRowEntryIndex = this->rowStarts[a_y+1];
		bool finished = false;

		for (int i=rowEntryIndex; i<nextRowEntryIndex && ~finished; ++i)
		{
			int columnNo = this->columnNos[i];
			std::cout << "cn" << columnNo << std::endl;

			// We've just gone past the column we want.
			if (columnNo > a_x)
			{
				typename std::vector<T>  ::iterator itMatrixEntries = this->matrixEntries.begin() + i-1;
				         std::vector<int>::iterator itColumnNos     = this->columnNos    .begin() + i-1;

				this->matrixEntries.insert(itMatrixEntries, a_value); // Is insert the best way of doing what I want? This will increase the lengths...
				this->columnNos    .insert(itColumnNos    , a_x);

				for (int j=a_y+1; j<=this->rowStarts.size(); ++j)
					++this->rowStarts[j];

				finished = true;
			}
		}*/
	}

	std::cout << "ROW STARTS" << std::endl;
	for (int i=0; i<this->rowStarts.size(); ++i)
		std::cout << rowStarts[i] << " ";
	std::cout << std::endl << std::endl;

	std::cout << "COLUMN NOS" << std::endl;
	for (int i=0; i<this->columnNos.size(); ++i)
		std::cout << columnNos[i] << " ";
	std::cout << std::endl << std::endl;

	std::cout << "MATRIX ENTRIES" << std::endl;
	for (int i=0; i<this->matrixEntries.size(); ++i)
		std::cout << matrixEntries[i] << " ";
	std::cout << std::endl << std::endl;
}


#endif