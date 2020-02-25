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
		// Resizing.
		virtual void resize(const int &a_noNonZeros) = 0; // Missing resize for number of rows and columns.

		// Gets an individual item.
		virtual const T item(const int &a_x, const int &a_y) const = 0;

	public:
		// Getters.
		virtual int    get_noRows() const = 0;
		virtual int    get_noColumns() const = 0;
		std::vector<T> get_diagonal() const;

		// Indexing.
		const T operator()(const int &a_x, const int &a_y) const;

		// Setters.
		virtual void set(const int &a_x, const int &a_y, const T &a_value) = 0;

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