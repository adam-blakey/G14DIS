/******************************************************************************
 * @details Declarations for [Matrix].
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/01/02
 ******************************************************************************/
#ifndef CLASS_MATRIX
#define CLASS_MATRIX

#include <cassert>
#include <vector>

template<class T>
class Matrix
{
	protected:
		// Storage.
		int noColumns;
		int noRows;

		// Resizing.
		virtual void resize(const int &a_noNonZeros) = 0; // Missing resize for number of rows and columns.

		// Gets an individual item.
		virtual T&       item(const int &a_x, const int &a_y) = 0;
		virtual const T& item(const int &a_x, const int &a_y) const = 0;

	public:
		// Getters.
		int            get_noRows() const;
		int            get_noColumns() const;
		std::vector<T> get_diagonal() const;

		// Indexing.
		T&       operator()(const int &a_x, const int &a_y);
		const T& operator()(const int &a_x, const int &a_y) const;

		// Matrix-Matrix operations.
		Matrix<T>& operator= (const Matrix<T> &a_RHS);
		Matrix<T>& operator+=(const Matrix<T> &a_RHS);
		Matrix<T>& operator-=(const Matrix<T> &a_RHS);
		Matrix<T>& operator*=(const Matrix<T> &a_RHS);

		// Matrix-scalar operations.
		Matrix<T>& operator*=(const T &a_RHS);
		Matrix<T>& operator/=(const T &a_RHS);

		// Matrix-vector operations.
		std::vector<T> operator*(const std::vector<T> &a_RHS) const;
};

#include "matrix.cpp"

#endif