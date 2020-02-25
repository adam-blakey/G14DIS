/******************************************************************************
 * @details Declarations for [Matrix].
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/01/02
 ******************************************************************************/
#ifndef CLASS_MATRIX_FULL
#define CLASS_MATRIX_FULL

#include "matrix.hpp"
#include <cassert>
#include <vector>

template<class T>
class Matrix_full: public Matrix<T>
{
	protected:
		// Hidden default constructor.
		Matrix_full();

		// Storage.
		std::vector<T> items;
		int noColumns;
		int noRows;

		// Resizing.
		void resize(const int &a_noNonZeros);

		// Index calculation.
		int get_index(const int &x, const int &y) const;

		// Gets an individual item.
			  T& item(const int &a_x, const int &a_y);
		const T  item(const int &a_x, const int &a_y) const;

	public:
		// Constructors.
		Matrix_full(const int &a_N);
		Matrix_full(const int &a_noColumns, const int &a_noRows);
		Matrix_full(const int &a_noColumns, const int &a_noRows, const T &a_initial);
		Matrix_full(const Matrix<T> &a_matrix);

		// Matrix-Matrix operations.
		Matrix_full<T> operator+(const Matrix<T> &a_RHS);
		Matrix_full<T> operator-(const Matrix<T> &a_RHS);
		Matrix_full<T> operator*(const Matrix<T> &a_RHS);

		// Matrix-scalar operations.
		Matrix_full<T> operator*(const T &a_RHS);
		Matrix_full<T> operator/(const T &a_RHS);

		// Getters.
		int get_noRows() const;
		int get_noColumns() const;
		
		// Setters.
		void set(const int &a_x, const int &a_y, const T &a_value);
};

#include "matrix_full.cpp"

#endif