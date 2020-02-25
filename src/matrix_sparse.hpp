/******************************************************************************
 * @details Declarations for [Matrix].
 * 
 * @author     Adam Matthew Blakey
 * @date       2020/01/02
 ******************************************************************************/
#ifndef CLASS_MATRIX_SPARSE
#define CLASS_MATRIX_SPARSE

#include "matrix.hpp"
#include <cassert>
#include <vector>

template<class T>
class Matrix_sparse: public Matrix<T>
{
	protected:
		// Hidden default constructor.
		Matrix_sparse();

		// Storage.
		std::vector<T>   matrixEntries;
		std::vector<int> columnNos;
		std::vector<int> rowStarts;
		int              noColumns;

		// Resizing.
		void resize(const int &a_noNonZeros);
		void resize(const int &a_noNonZeros, const int &a_noRows);
		void resize(const int &a_noNonZeros, const int &a_noRows, const int &a_noColumns);

		// Index calculation.
		int get_index(const int &a_x, const int &a_y) const;

		// Gets an individual item.
		const T item(const int &a_x, const int &a_y) const;

	public:
		// Constructors.
		Matrix_sparse(const int &a_noNonZero, const int &a_noRows, const int &a_noColumns);
		Matrix_sparse(const Matrix<T> &a_matrix);
		//add constructor with columnNos and rowStarts already added?

		// Matrix-Matrix operations.
		Matrix_sparse<T> operator+(const Matrix<T> &a_RHS);
		Matrix_sparse<T> operator-(const Matrix<T> &a_RHS);
		Matrix_sparse<T> operator*(const Matrix<T> &a_RHS);

		// Matrix-scalar operations.
		Matrix_sparse<T> operator*(const T &a_RHS);
		Matrix_sparse<T> operator/(const T &a_RHS);

		// Getters.
		int get_noNonZero() const;
		int get_noRows() const;
		int get_noColumns() const;

		// Setting.
		void set(const int &a_x, const int &a_y, const T &a_value);
};

#include "matrix_sparse.cpp"

#endif